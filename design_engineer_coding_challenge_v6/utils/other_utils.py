from io import StringIO
import csv

import pysam
from Bio import SeqIO


def write_variant_collection_to_file(
    variant_collection, output_file, collection_name, species
):
    with open(output_file, "w") as csv_file:
        header_line_1 = [
            "DesignCollectionName",
            collection_name,
            *[None] * 9,
        ]  # 12 columns in total
        header_line_2 = [
            "CodonUsage",
            species,
            *[None] * 9,
        ]  # 12 columns in total
        header_line_3 = [
            "EditType",
            "TargetType",
            "FeatureType",
            "TargetName",
            "CoordinateType",
            "PositionType",
            "PositionValue",
            "SaturationReplacementType",
            "Control",
            "NumberOfCoordinatesToDelete",
            "InsertionSequence",
            "Comment",
        ]

        writer = csv.writer(csv_file)
        writer.writerows([header_line_1, header_line_2, header_line_3])

        # TODO Generalize this csv writing part; avoid hard coding?
        TargetType = "Chromosome"
        FeatureType = ""
        CoordinateType = "Nucleotide"
        PositionType = "Specific Position"
        SaturationReplacementType = ""
        Control = ""
        Comment = ""
        NumberOfCoordinatesToDelete = ""
        InsertionSequence = ""
        for variants in variant_collection:
            TargetName = f"chr{variants.contig}"
            PositionValue = variants.start
            if variants.is_deletion:
                EditType = "Deletion"
                NumberOfCoordinatesToDelete = len(variants.ref)
            elif variants.is_insertion:
                EditType = "Insertion"
                InsertionSequence = variants.alt
            else:
                EditType = "Substitution"
                InsertionSequence = variants.alt

            writer.writerow(
                [
                    EditType,
                    TargetType,
                    FeatureType,
                    TargetName,
                    CoordinateType,
                    PositionType,
                    PositionValue,
                    SaturationReplacementType,
                    Control,
                    NumberOfCoordinatesToDelete,
                    InsertionSequence,
                    Comment,
                ]
            )


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

    # logging.info(
    #     f"Retrieving sequence from {fasta_file_path}, region {sequence_header_region}"
    # )

    try:
        fasta_sequence_record = pysam.faidx(fasta_file_path, sequence_header_region)
        sequence_record = SeqIO.read(StringIO(fasta_sequence_record), "fasta").seq
    except Exception as e:
        raise f"Failed to retrieve sequence from FASTA file {fasta_file_path} with error {e}"
    return sequence_record
