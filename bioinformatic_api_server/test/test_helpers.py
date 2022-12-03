import pytest
from service.utility_service import reverse_complement, parse_genome_data
from service.query_handler_service import (
    get_genomes,
    get_length,
    searchseq,
    retrieveseq,
    read_index_file,
)

from service.genome_handler_sevice import register_genome


@pytest.mark.parametrize(
    "sequence, exp_revcomped",
    [
        ("aaccttgca", "tgcaaggtt"),
        ("AAccttgca", "tgcaaggTT"),
        ("AAccNtgca", "tgcaNggTT"),
    ],
)
def test_reverse_complement(sequence, exp_revcomped):
    """F"""
    obs_revcomped = reverse_complement(sequence)
    assert obs_revcomped == exp_revcomped


def test_parse_genome_data():
    contigs_and_seqs = {
        "seq_1": "ATGAAAGCGAACCTGCTGGTTCTGCTGTGCGCGCTGGCGGCGGCGGACGCGGACACCATCTGCATCGGTTACCACGCGAACAACTCTACCGACACCGTTGACACCGTTCTGGAAAAAAACGTTACCGTTACCCACTCTGTTAACCTGCTGGAAGACTCTCACAACGGTAAACTGTGCCGT",
        "seq_2": "CTGAAAGGTATCGCGCCGCTGCAGCTGGGTAAATGCAACATCGCGGGTTGGCTGCTGGGT",
        "seq_3": "AACCCGGAATGCGACCCGCTGCTGCCGGTTCGTTCTTGGTCTTACATCGTTGAAACCCCGAACTCTGAAAACGGTATCTGCTACCCGGGTGACTTCATCGACTACGAAGAACTGCGTGAA",
    }

    def create_genome_data(contigs_and_seqs, length_seq_line=80):
        genome_data = ""
        for contig, seq in contigs_and_seqs.items():
            genome_data += f"> {contig}\n"
            genome_data += (
                "\n".join(
                    [
                        seq[idx : idx + length_seq_line]
                        for idx in range(0, len(seq), length_seq_line)
                    ]
                )
                + "\n"
            )
        return genome_data

    genome_data = create_genome_data(contigs_and_seqs)
    obs_contigs_and_seqs = parse_genome_data(genome_data)
    assert contigs_and_seqs == obs_contigs_and_seqs


@pytest.fixture()
def fasta_file_path():
    return "test/data/test_fasta.fa"


def test_get_genomes():
    genome_register = "test/data/genome_register.csv"
    genome_register_data = {
        "18967147e58946829547cfb9a897af4f": {
            "unique_identifier": "18967147e58946829547cfb9a897af4f",
            "upload_path": "data/18967147e58946829547cfb9a897af4f.fa",
            "upload_name": "sample_1.fasta",
        }
    }
    assert genome_register_data == get_genomes(genome_register)


@pytest.mark.parametrize(
    "sequence_header_region, expected_length",
    [
        ("", {"test_sequence": "368"}),
        ("test_sequence", {"test_sequence": "368"}),
    ],
)
def test_get_length(sequence_header_region, expected_length, fasta_file_path):
    assert expected_length == get_length(fasta_file_path, sequence_header_region)


def test_read_index_file(fasta_file_path):
    index_data = {"test_sequence": "368"}
    assert index_data == read_index_file(fasta_file_path)


def test_retrieveseq(fasta_file_path):
    expected_sequence = ">test_sequence:1-10\nACAAGATGCC\n"
    sequence_header = "test_sequence:1-10"
    assert expected_sequence == retrieveseq(fasta_file_path, sequence_header)


def test_searchseq(fasta_file_path):
    search_sequence = "TCC"
    sequence_header = "test_sequence"
    expected_position = {
        "test_sequence": {
            "forward_direction": "15-17",
            "reverse_compliment_direction": "69-71",
        }
    }
    assert expected_position == searchseq(
        fasta_file_path, sequence_header, search_sequence
    )
