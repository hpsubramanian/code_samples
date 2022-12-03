import os

from handlers.base_handler import BaseView
from service import query_handler_service
from config import upload_folder


class QueryEngine(BaseView):
    """Query genomes."""

    SUPPORTED_METHODS = "GET"

    def get(self, query_type):
        """Routes requests based on query_type parsed from url"""

        unique_identifier = self.get_query_argument("uid", None)
        sequence_header_region = self.get_query_argument("sequence_region", "")
        query_sequence = self.get_query_argument("query", None)

        fasta_file_name = f"{unique_identifier}.fa"
        fasta_file_path = os.path.join(upload_folder, fasta_file_name)

        if query_type == "listgenomes":
            genome_data = query_handler_service.get_genomes()
            if genome_data:
                self.send_response(genome_data)
            return

        if (
            query_type in ["length", "retrieveseq", "searchseq"]
            and unique_identifier is None
        ):
            raise ValueError("unique_identifier is a required parameter")

        if query_type == "length":
            try:
                unique_identifier = unique_identifier.strip("/")
                sequence_header_region = sequence_header_region.strip("/")
                length_data = query_handler_service.get_length(
                    fasta_file_path, sequence_header_region
                )
                if length_data:
                    self.send_response(length_data)
                else:
                    self.send_response("No length data was returned")
            except Exception as e:
                raise f"Failed to get sequence length - {e}"

            return

        if query_type == "retrieveseq":
            try:
                sequence_info = query_handler_service.retrieveseq(
                    fasta_file_path, sequence_header_region
                )
                if sequence_info:
                    self.send_response(sequence_info)
                return
            except Exception as e:
                raise f"Failed to get sequence information - {e}"

        if query_type == "searchseq":
            try:
                searchseq_info = query_handler_service.searchseq(
                    fasta_file_path, sequence_header_region, query_sequence
                )
                self.send_response(searchseq_info)
                return
            except Exception as e:
                raise f"Failed to get search for sequence - {e}"
