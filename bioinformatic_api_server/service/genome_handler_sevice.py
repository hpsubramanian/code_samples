import uuid
from collections import defaultdict
import pysam
import os.path
import logging

from config import genome_register, genome_register_headers, upload_folder, log_level

logger = logging.getLogger(__name__)
logger.setLevel(log_level)


def register_genome(request_file_items):
    """
    Function to store, index and register genomes
    Parameters
    ----------
    request_file_items
        Files along with metadata from request

    Returns
    -------
    Dictionary of lists with unique identifier and fasta metadata

    """
    logging.info("Received request to register fasta files")
    # TODO FASTA File extension check, file content/format check
    # TODO Separate FASTA index creation - useful for large FASTA files
    # TODO Have a better way to handle genome register - some DB, more metadata, md5 hash
    # TODO Auto create upload folders during startup

    if not os.path.isfile(genome_register):
        logging.info("Writing genome register headers")
        write_content_to_file(genome_register, ",".join(genome_register_headers))
        write_content_to_file(genome_register, "\n")

    # Adapted from Tornado docs
    uploaded_files = defaultdict(list)
    for field_names, files in request_file_items:
        for info in files:
            filename = info["filename"]
            body = info["body"].decode("utf8")

            logging.info(f"Processing {filename}")

            # Determine file name, upload path
            unique_filename = uuid.uuid4().hex
            upload_path = f"{upload_folder}/{unique_filename}.fa"
            genome_register_entry = f"{unique_filename},{upload_path},{filename}\n"

            write_content_to_file(file_path=upload_path, content=body)
            write_content_to_file(
                file_path=genome_register, content=genome_register_entry
            )

            pysam.faidx(upload_path)
            uploaded_files["uploaded_files"].append(
                {"filename": filename, "unique_identifier": unique_filename}
            )
    return uploaded_files


def write_content_to_file(file_path, content, mode="a+"):
    """
    Function to write contents to a file
    Parameters
    ----------
    file_path
        file to write data to
    content
        data to write to file_path
    mode
        mode to open the file in

    Returns
    -------

    """
    logging.info(f"Writing content to {file_path}")
    try:
        with open(file_path, mode) as in_file:
            in_file.write(content)
    except Exception as e:
        raise f"Error writing content to file- {e}"
