from handlers.base_handler import BaseView
from service import genome_handler_sevice


class GenomeHandler(BaseView):
    """Register genomes."""

    SUPPORTED_METHODS = ("POST",)

    def post(self):
        """Receives FASTA files to register them"""
        unique_filenames = genome_handler_sevice.register_genome(
            self.request.files.items()
        )
        self.send_response(unique_filenames)
