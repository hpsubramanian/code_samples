import logging
from itertools import product

from pyensembl import ensembl_grch38
from varcode import Variant


def return_bases_of_length_l(length):
    """Function to return all nucleotide combinations/kmers of length l
    Adapted from
    https://stackoverflow.com/questions/48677692/generating-all-possible-k-mers-string-combinations-from-a-given-list
    Parameters
    ----------
    length : int
        length of the nucleotide combinations

    Returns
    -------
    list of all possible nucleotide combinations/kmers of a given length
    """

    return ["".join(base) for base in product("ATGC", repeat=length)]


def genomic_positions_list(exon_start, exon_end):
    """
    Function to return genomic positions as a list, given a start and stop
    Parameters
    ----------
    exon_start : int
        exon start position
    exon_end : int
        exon end position

    Returns
    -------
    genomic_positions_in_exon : list
        list of all genomic positions between start and stop

    """
    genomic_positions_in_exon = [
        number for number in range(exon_start, exon_end + 1)
    ]  # Add +1 to end to include the end base
    return genomic_positions_in_exon


def generate_substitutions(
    chromosome,
    exon_start,
    exon_end,
    exon_sequence,
    substitution_length,
    allowed_effects,
):
    """
    Function to generate all possible substitutions of given length in an exon
    Parameters
    ----------
    chromosome : int
        chromosome number of the exon
    exon_start : int
        exon start position
    exon_end : int
        exon stop position
    exon_sequence : str
        exon nucleotide sequence
    substitution_length : int
        length of substitution variant
    allowed_effects : list
        list of allowed variant effects, used as a filter

    Returns
    -------
    substitution_variants : list
        a list of substitution variants.
    """
    substitution_variants = []
    genomic_positions_in_exon = genomic_positions_list(
        exon_start=exon_start, exon_end=exon_end
    )
    alternate_kmers = return_bases_of_length_l(length=substitution_length)
    index = 0
    for genomic_positions in genomic_positions_in_exon:
        ref_base = exon_sequence[index : index + substitution_length]
        for alternate_base in alternate_kmers:
            current_variant = Variant(
                contig=chromosome,
                start=genomic_positions,
                ref=ref_base,
                alt=alternate_base,
                ensembl=ensembl_grch38,
            )
            if (
                type(current_variant.effects().top_priority_effect()).__name__
                in allowed_effects
            ):
                substitution_variants.append(current_variant)
        # TODO Possible Logic issue here - need to check how it handles last bases or intronic sequences
        if index + substitution_length >= len(exon_sequence):
            break
        index += 1
    return substitution_variants


def generate_insertions(
    chromosome, exon_start, exon_end, exon_sequence, insertion_length, allowed_effects
):
    """
    Function to generate all possible insertions of given length in an exon
    Parameters
    ----------
    chromosome : int
        chromosome number of the exon
    exon_start : int
        exon start position
    exon_end : int
        exon stop position
    exon_sequence : str
        exon nucleotide sequence
    substitution_length : int
        length of substitution variant
    allowed_effects : list
        list of allowed variant effects, used as a filter

    Returns
    -------
    insertion_variants : list
        a list of insertion variants.

    """
    insertion_variants = []
    genomic_positions_in_exon = genomic_positions_list(
        exon_start=exon_start, exon_end=exon_end
    )
    alternate_kmers = return_bases_of_length_l(length=insertion_length)
    for position, ref_base in zip(genomic_positions_in_exon, exon_sequence):
        for alternate_base in alternate_kmers:
            current_variant = Variant(
                contig=chromosome,
                start=position,
                ref=ref_base,
                alt=ref_base + alternate_base,
                ensembl=ensembl_grch38,
            )
            if (
                type(current_variant.effects().top_priority_effect()).__name__
                in allowed_effects
            ):
                insertion_variants.append(current_variant)
    return insertion_variants


def generate_deletions(
    chromosome, exon_start, exon_end, exon_sequence, deletion_length, allowed_effects
):
    """
    Function to generate all possible deletions of given length in an exon
    Parameters
    ----------
    chromosome : int
        chromosome number of the exon
    exon_start : int
        exon start position
    exon_end : int
        exon stop position
    exon_sequence : str
        exon nucleotide sequence
    substitution_length : int
        length of substitution variant
    allowed_effects : list
        list of allowed variant effects, used as a filter

    Returns
    -------
    deletion_variants : list
        a list of deletion variants.
    """
    deletion_variants = []
    index = 0
    genomic_positions_in_exon = genomic_positions_list(
        exon_start=exon_start, exon_end=exon_end
    )
    for genomic_positions in genomic_positions_in_exon:
        ref_base = exon_sequence[
            index : index + deletion_length + 1
        ]  # Deletion, hence we want to get an additional base
        alt_base = exon_sequence[index]
        current_variant = Variant(
            contig=chromosome,
            start=genomic_positions,
            ref=ref_base,
            alt=alt_base,
            ensembl=ensembl_grch38,
        )
        # TODO Take this out of this function; find a better way instead of iterating twice if taken out
        if (
            type(current_variant.effects().top_priority_effect()).__name__
            in allowed_effects
        ):
            deletion_variants.append(current_variant)
        if index + deletion_length + 1 > len(exon_sequence):
            break
        index += 1
    return deletion_variants
