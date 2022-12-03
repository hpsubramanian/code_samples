from utils import variant_utils


def test_return_bases_of_length_l():
    length_one = ["A", "T", "G", "C"]
    length_two = [
        "AA",
        "AT",
        "AG",
        "AC",
        "TA",
        "TT",
        "TG",
        "TC",
        "GA",
        "GT",
        "GG",
        "GC",
        "CA",
        "CT",
        "CG",
        "CC",
    ]
    length_three = [
        "AAA",
        "AAT",
        "AAG",
        "AAC",
        "ATA",
        "ATT",
        "ATG",
        "ATC",
        "AGA",
        "AGT",
        "AGG",
        "AGC",
        "ACA",
        "ACT",
        "ACG",
        "ACC",
        "TAA",
        "TAT",
        "TAG",
        "TAC",
        "TTA",
        "TTT",
        "TTG",
        "TTC",
        "TGA",
        "TGT",
        "TGG",
        "TGC",
        "TCA",
        "TCT",
        "TCG",
        "TCC",
        "GAA",
        "GAT",
        "GAG",
        "GAC",
        "GTA",
        "GTT",
        "GTG",
        "GTC",
        "GGA",
        "GGT",
        "GGG",
        "GGC",
        "GCA",
        "GCT",
        "GCG",
        "GCC",
        "CAA",
        "CAT",
        "CAG",
        "CAC",
        "CTA",
        "CTT",
        "CTG",
        "CTC",
        "CGA",
        "CGT",
        "CGG",
        "CGC",
        "CCA",
        "CCT",
        "CCG",
        "CCC",
    ]

    assert sorted(length_one) == sorted(variant_utils.return_bases_of_length_l(1))
    assert sorted(length_two) == sorted(variant_utils.return_bases_of_length_l(2))
    assert sorted(length_three) == sorted(variant_utils.return_bases_of_length_l(3))


def test_genomic_positions_list():
    start = 10
    end = 20
    positions = [x for x in range(start, end + 1)]
    assert positions == variant_utils.genomic_positions_list(10, 20)


def test_generate_substitutions():
    all_possible_effects = [
        "AlternateStartCodon",
        "ComplexSubstitution",
        "Deletion",
        "ExonLoss",
        "ExonicSpliceSite",
        "FivePrimeUTR",
        "FrameShift",
        "FrameShiftTruncation",
        "IncompleteTranscript",
        "Insertion",
        "Intergenic",
        "Intragenic",
        "Intronic",
        "IntronicSpliceSite",
        "NoncodingTranscript",
        "PrematureStop",
        "Silent",
        "SpliceAcceptor",
        "SpliceDonor",
        "StartLoss",
        "StopLoss",
        "Substitution",
        "ThreePrimeUTR",
    ]
    chromosome = 17
    exon_start = 43063874
    exon_end = 43063951
    exon_sequence = (
        "AGAAATAGCTAACTACCCATTTTCCTCCCGCAATTCCTAGAAAATATTTCAGTGTCCGTTCACACACAAACTCAGCAT"
    )
    total_number_of_possible_1b_substitutions = (
        len(exon_sequence) * 4
    )  # Each position can have all four bases;
    assert total_number_of_possible_1b_substitutions == len(
        variant_utils.generate_substitutions(
            chromosome=chromosome,
            exon_start=exon_start,
            exon_end=exon_end,
            allowed_effects=all_possible_effects,
            substitution_length=1,
            exon_sequence=exon_sequence,
        )
    )
