import re
import csv
import pysam
import logging
from io import StringIO

from Bio import SeqIO
from config import genome_register, log_level

logger = logging.getLogger(__name__)
logger.setLevel(log_level)


def get_genomes(genome_file=genome_register):
    """
    Function to return all available/registered genomes

    Parameters
    ----------
    genome_file
        Genome register containing all available genomes

    Returns
    -------
    dict
        A dict of unique genome identifiers, upload path and upload name.

    """
    logging.info("Fetching genomes from register")

    data = {}
    with open(genome_file, "r") as genome_register_fh:
        genome_entries = csv.DictReader(genome_register_fh)
        for rows in genome_entries:
            key = rows["unique_identifier"]
            data[key] = rows

    logging.info("Returning genomes available in the genome register")

    return data


def get_length(fasta_file_path, sequence_header_region):
    """
    Function to get length of sequence(s) in a fasta file

    Parameters
    ----------
    fasta_file_path
        Path to fasta file whose length is being queried
    sequence_header_region
        sequence name whose length is required

    Returns
    -------
    length_dictionary
        Dictionary containing sequence names and lengths

    """

    logging.info(
        f"Fetching length of FASTA Sequences: Identifier {fasta_file_path} Region {sequence_header_region}"
    )

    length_dictionary = read_index_file(fasta_file_path)
    if sequence_header_region:
        try:
            chromosome_length_info = length_dictionary[sequence_header_region]
            return {sequence_header_region: chromosome_length_info}
        except Exception as e:
            raise f"Error finding length of {sequence_header_region} sequence - {e}"
    else:
        return length_dictionary


def read_index_file(fasta_file_path):
    """
    Function to read FASTA index file

    Parameters
    ----------
    fasta_file_path
        Path of fasta file whose index needs to be read

    Returns
    -------
    length_dictionary
        Dictionary containing sequence names and lengths

    """

    logging.info(f"Reading FASTA index file for {fasta_file_path}")

    index_file = f"{fasta_file_path}.fai"
    length_dictionary = {}
    with open(index_file, "r") as index_fh:
        for line in index_fh:
            sequence_name, length, *unwanted = line.strip("\n").split("\t")
            length_dictionary[sequence_name] = length
    return length_dictionary


def retrieveseq(fasta_file_path, sequence_header_region):
    """
    Get subsequence from a fasta file

    Parameters
    ----------
    fasta_file_path
        FASTA file path to query
    sequence_header_region
        Sequence_name:start-stop, can be Sequence_name alone as well

    Returns
    -------
    sequence_record
        Sequences at the requested positions

    """

    logging.info(
        f"Retrieving sequence from {fasta_file_path}, region {sequence_header_region}"
    )

    try:
        sequence_record = pysam.faidx(fasta_file_path, sequence_header_region)
    except Exception as e:
        raise f"Failed to retrieve sequence from FASTA file {fasta_file_path} with error {e}"
    return sequence_record


# TODO handle with base exception
def searchseq(fasta_file_path, sequence_header_region, search_sequence):
    """
    Searches for a substring in FASTA file and returns coordinates.

    Searches in two directions - forward and reverse complement

    Parameters
    ----------
    fasta_file_path
         FASTA file path to search
    sequence_header_region
        Sequence_name:start-stop, can be Sequence_name alone as well
    search_sequence
        Query sequence

    Returns
    -------
    sequence_location
        Dictionary containing matches in forward and reverse complement direction

    """

    logging.info(
        f"Searching for sequence {search_sequence} in {fasta_file_path}, region {sequence_header_region}"
    )

    query_sequence = retrieveseq(fasta_file_path, sequence_header_region)

    # Convert fasta to a SeqIO record and get raw sequence
    fasta_sequence_record = SeqIO.read(StringIO(query_sequence), "fasta").seq
    length_of_search_sequence = len(search_sequence)

    # Get start and end position from sequence_header_region
    regex = r"(\w*)(:(\d*-\d*))?"
    matches = re.search(regex, sequence_header_region)
    start_stop = matches.group(3)
    if start_stop:
        start, end = start_stop.split("-")
    else:
        start = 1

    # Find position of sequence in query in both directions
    start_in_forward_dir, start_in_reverse_com_dir = search_seq_both_dir(
        fasta_sequence_record, search_sequence
    )

    actual_fwd_start, actual_fwd_end = None, None
    actual_rc_start, actual_rc_end = None, None

    if start_in_forward_dir != -1:
        actual_fwd_start = start_in_forward_dir + int(start)
        actual_fwd_end = actual_fwd_start + length_of_search_sequence - 1

    if start_in_reverse_com_dir != -1:
        actual_rc_start = start_in_reverse_com_dir + int(start)
        actual_rc_end = actual_rc_start + length_of_search_sequence - 1

    sequence_location = {
        sequence_header_region: {
            "forward_direction": f"{actual_fwd_start}-{actual_fwd_end}",
            "reverse_compliment_direction": f"{actual_rc_start}-{actual_rc_end}",
        }
    }

    return sequence_location


def search_seq_both_dir(fasta_sequence_record, search_sequence):
    """
    Find start position of a substring in a fasta file

    Parameters
    ----------
    fasta_sequence_record
        Target Sequence record to search
    search_sequence
        Query sequence

    Returns
    -------
    start_in_forward_dir
        Position of match in forward direction
    start_in_reverse_com_dir
        Position of match in reverse complement direction

    """

    logging.info(f"Searching for sequence {search_sequence}")

    start_in_forward_dir = fasta_sequence_record.find(search_sequence)
    start_in_reverse_com_dir = fasta_sequence_record.reverse_complement().find(
        search_sequence
    )

    return start_in_forward_dir, start_in_reverse_com_dir
