import click
import logging
import os

from pyensembl import EnsemblRelease
from utils import variant_utils, other_utils, transcript_utils
from varcode import VariantCollection

logging.getLogger("varcode").setLevel(logging.WARNING)
logging.getLogger("pyensembl").setLevel(logging.WARNING)
os.environ["PYENSEMBL_CACHE_DIR"] = "./data"


def generate_variants(
    gene_name,
    variation_length_list,
    allowed_effects,
    substitutions_flag,
    insertions_flag,
    deletions_flag,
):
    """Function to generate variant in a given gene, of given length and type

    Parameters
    ----------
    gene_name : str
        HGNC approve gene name to generate variants in
    variation_length_list : list
        list containing lengths of variations to generate
    allowed_effects : list
        list of allowed variant effects; see varcode docs for allowed terms
    substitutions_flag : bool
        Boolean flag to generate substitutions
    insertions_flag : bool
        Boolean flag to generate insertions
    deletions_flag : bool
        Boolean flag to generate deletions

    Returns
    -------
    VariantCollection
        Variant Collection object containing all variants generated of desired length and type
    """

    ensembl_handler = EnsemblRelease(release=104, species="Human")
    logging.info("Using EnsemblRelease 104, species human")

    logging.info(f"Generating variants in {gene_name}")
    gene_object_list = ensembl_handler.genes_by_name(gene_name)

    # Ensure only one gene object is returned by ensembl_handler
    if len(gene_object_list) != 1:
        raise ValueError(
            "More than 1 gene was found; Provide ensembl compatible, HGNC approved gene name"
        )
    else:
        gene_object = gene_object_list[0]
        gene_id = gene_object.gene_id

    # Find ensembl transcript id and exons in canonical transcript
    canonical_transcript = transcript_utils.get_canonical_transcript(
        ensembl_gene_id=gene_id
    )
    logging.info(f"The canonical transcript for {gene_name} is {canonical_transcript}")
    canonical_transcript_exons = transcript_utils.get_exons_for_transcripts(
        pyensembl_handler=ensembl_handler, transcript_id_list=[canonical_transcript]
    )

    # Find ensembl transcript id and exons in all protein coding transcripts
    complete_protein_coding_transcripts = (
        transcript_utils.get_complete_protein_coding_transcripts(
            gene_object=gene_object
        )
    )
    exons_in_all_protein_coding_transcripts = (
        transcript_utils.get_exons_for_transcripts(
            pyensembl_handler=ensembl_handler,
            transcript_id_list=complete_protein_coding_transcripts,
        )
    )

    # Find exon usage across different transcripts
    exon_usage_counter = transcript_utils.count_exon_usage(
        exon_id_list=exons_in_all_protein_coding_transcripts
    )

    # Take top 20% of the most used transcripts
    # TODO Rename percentage_value
    # TODO Get value from command line
    most_commonly_used_exons = transcript_utils.find_top_n_exons(
        exon_counter=exon_usage_counter, percentage_value=0.2
    )

    # Combine both exon lists to be comprehensive
    list_of_exons_to_introduce_variants = list(
        set(canonical_transcript_exons + most_commonly_used_exons)
    )

    master_variant_list = []
    exon_index = 0

    # Iterate through exons, iterate through variation length, generate a particular type if the flag is set
    # Filter by allowed effects
    for exon_identifier in list_of_exons_to_introduce_variants:
        exon_index += 1
        logging.info(
            f"Processing exon {exon_index} of {len(list_of_exons_to_introduce_variants)}"
        )

        # Get exon information
        exon_start = ensembl_handler.exon_by_id(exon_identifier).start
        exon_end = ensembl_handler.exon_by_id(exon_identifier).end
        chromosome = ensembl_handler.exon_by_id(exon_identifier).contig
        fasta_file_path = f"./data/Homo_sapiens.GRCh38.dna.chromosome.{chromosome}.fa"
        fasta_sequence_record = str(
            other_utils.retrieveseq(
                fasta_file_path, f"{chromosome}:{exon_start}-{exon_end}"
            )
        )
        for variation_length in variation_length_list:

            substitutions, insertions, deletions = [], [], []
            if substitutions_flag:
                logging.debug(
                    f"Generating {variation_length}bp substitutions in {exon_identifier}"
                )
                substitutions = variant_utils.generate_substitutions(
                    chromosome=chromosome,
                    exon_start=exon_start,
                    exon_end=exon_end,
                    exon_sequence=fasta_sequence_record,
                    substitution_length=variation_length,
                    allowed_effects=allowed_effects,
                )
                logging.debug(f"In {exon_identifier}, generated {len(substitutions)}")
                master_variant_list.extend(substitutions)

            if insertions_flag:
                logging.debug(
                    f"Generating {variation_length}bp insertions in {exon_identifier}"
                )
                insertions = variant_utils.generate_insertions(
                    chromosome=chromosome,
                    exon_start=exon_start,
                    exon_end=exon_end,
                    exon_sequence=fasta_sequence_record,
                    insertion_length=variation_length,
                    allowed_effects=allowed_effects,
                )
                logging.debug(f"In {exon_identifier}, generated {len(insertions)}")
                master_variant_list.extend(insertions)

            if deletions_flag:
                logging.debug(
                    f"Generating {variation_length}bp deletions in {exon_identifier}"
                )
                deletions = variant_utils.generate_deletions(
                    chromosome=chromosome,
                    exon_start=exon_start,
                    exon_end=exon_end,
                    exon_sequence=fasta_sequence_record,
                    deletion_length=variation_length,
                    allowed_effects=allowed_effects,
                )
                logging.debug(f"In {exon_identifier}, generated {len(deletions)}")
                master_variant_list.extend(deletions)

    return VariantCollection(variants=master_variant_list)


@click.command()
@click.option(
    "--target_type",
    required=True,
    type=click.Choice(["gene", "transcript", "exon"], case_sensitive=False),
    help="Type of annotation target to introduce variants in",
)
@click.option(
    "--gene_name",
    required=True,
    type=str,
    help="HGNC Approved Gene Name. Required if target_type is gene, ignored otherwise",
)
@click.option(
    "--variation_length",
    default="1",
    help="Lengths of Insertions,Deletions,Substitutions to simulate. Comma seperated list",
    show_default=True,
)
@click.option(
    "--substitutions/--no-substitutions",
    default=True,
    help="Flag to generate substitutions",
    show_default=True,
)
@click.option(
    "--insertions/--no-insertions",
    default=False,
    help="Flag to generate insertions",
    show_default=True,
)
@click.option(
    "--deletions/--no-deletions",
    default=False,
    help="Flag to generate deletions",
    show_default=True,
)
@click.option(
    "--output_file_path", type=str, help="Path for output design file to write"
)
def main(
    target_type,
    gene_name,
    variation_length,
    substitutions,
    insertions,
    deletions,
    output_file_path,
):
    """
    Program to generate variants of a desired effect in a human gene/transcript/exon.
    Run variant_generator.py --help for information on input parameters.
    Default path to output design file is "Human_{gene_name}_variants.csv"
    """

    allowed_effects = ["ExonLoss", "FrameShiftTruncation", "PrematureStop"]
    variation_length_list = [int(length) for length in variation_length.split(",")]

    variant_collection = generate_variants(
        gene_name=gene_name,
        variation_length_list=variation_length_list,
        allowed_effects=allowed_effects,
        substitutions_flag=substitutions,
        insertions_flag=insertions,
        deletions_flag=deletions,
    )

    if output_file_path is None:
        output_file_path = f"Human_{gene_name}_variants.csv"

    other_utils.write_variant_collection_to_file(
        variant_collection=variant_collection,
        output_file=output_file_path,
        species="Human",
        collection_name=f"Human_{gene_name}_variants",
    )


if __name__ == "__main__":
    main()
