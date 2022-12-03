from Bio import SeqIO
from Bio.Seq import Seq
from typing import Dict
from io import StringIO


def reverse_complement(sequence: str) -> str:
    """Reverse complement a sequence."""
    revcomped = Seq(sequence).reverse_complement()
    return str(revcomped)


def parse_genome_data(contents: str) -> Dict[str, str]:
    """Parse genome data.

    Parameters
    ----------
    contents
        Contents of the genome file data.

    Returns
    -------
    Dictionary mapping names to sequences.

    """
    buffer = StringIO(contents)
    return {record.name: str(record.seq) for record in SeqIO.parse(buffer, "fasta")}
