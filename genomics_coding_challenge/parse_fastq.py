import logging
import click
import os

from utilities import plotting_utils, umi_utils
import pandas as pd


@click.command()
@click.option(
    "--input_fastq",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ file to extract UMI from. Can be gzipped",
)
@click.option(
    "--umi_prefix",
    required=True,
    type=str,
    help="UMI Prefix",
)
@click.option(
    "--umi_suffix",
    required=True,
    type=str,
    help="UMI Suffix",
)
@click.option(
    "--dendrogram_clustering_method",
    required=True,
    type=click.Choice(
        ["single", "complete", "average", "weighted"],
        case_sensitive=False,
    ),
    default="single",
    help="Method to be used for clustering while generating a dendrogram, based on Hamming distances. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage",
)
def main(input_fastq, umi_prefix, umi_suffix, dendrogram_clustering_method):

    # Make output directory if it doesn't exist
    if not os.path.exists("./data"):
        os.makedirs("./data")

    # Read input fastq file and extract UMI's. Write dataframe to file.
    # Note: Probably unnecessary to save sequence headers -> UMI information. Just UMI's and counts should do in this case
    main_readname_umi_df = umi_utils.extract_umi_from_fastq(
        input_fastq, umi_prefix, umi_suffix
    )
    umi_utils.write_data_frame_to_file(
        dataframe=main_readname_umi_df, filename="./data/UMI_in_each_read.txt"
    )

    # Get unique UMI counts and write dataframe to file.
    unique_umi_counts = umi_utils.get_unique_umi_counts(main_readname_umi_df)
    unique_umi_counts_df = pd.DataFrame.from_records(
        unique_umi_counts.most_common(), columns=["UMI_Sequence", "Count"]
    )
    umi_utils.write_data_frame_to_file(
        dataframe=unique_umi_counts_df, filename="./data/Unique_UMI_Sequence_Counts.txt"
    )

    logging.info(f"There were {len(unique_umi_counts)} unique UMI's in the FASTQ file")

    # Plot UMI Count Frequency distribution and write dataframe to file
    umi_count_frequency_table = plotting_utils.create_umi_frequency_distribution_graph(
        unique_umi_counts, "./data/UMI_Count_Frequency_distribution.png"
    )
    umi_utils.write_data_frame_to_file(
        dataframe=umi_count_frequency_table,
        filename="./data/UMI_Count_Frequency_Table.txt",
    )

    # Calculating pairwise hamming distance between UMI's and saving it to a file.
    (
        pdist_distance_matrix,
        squared_distance_matrix,
    ) = umi_utils.calculate_pairwise_hamming_distance(unique_umi_counts)

    umi_utils.write_data_frame_to_file(
        dataframe=squared_distance_matrix,
        filename="./data/UMI_Hamming_Distance_Table.txt",
        write_index=True,
    )

    # Obtain labels to be used in the dendrogram
    labels = [
        f"{umi} ({unique_umi_counts[umi]})" for umi in sorted(unique_umi_counts.keys())
    ]

    # Plot dendrogram
    plotting_utils.create_dendrogram(
        distance_matrix=pdist_distance_matrix,
        clustering_method=dendrogram_clustering_method,
        output_file_name=f"./data/UMI_Dendrogram_{dendrogram_clustering_method}.png",
        labels=labels,
    )

    # Find similar UMI's and calculate incorrect bases
    logging.info("Finding similar UMI's and calculating incorrectly sequenced bases")

    similar_umis = umi_utils.find_similar_umis(
        hamming_distance_matrix=squared_distance_matrix, distance_value=2
    )
    error_umi_bases = umi_utils.collapse_umis_find_error_bases(
        unique_umi_counts=unique_umi_counts, similar_umis_list_of_lists=similar_umis
    )

    # Find total UMI bases
    total_umi_bases = umi_utils.calculate_total_umi_bases_sequenced(unique_umi_counts)

    # Calculate and log % of bases called accurately.
    accurate_base_calling_percentage = (
        total_umi_bases - error_umi_bases
    ) / total_umi_bases

    logging.info(
        "Percentage of accurately called bases - (number of correctly called bases / total called bases)*100"
    )
    logging.info(
        f"Percentage of accurately called bases - ({total_umi_bases - error_umi_bases} / {total_umi_bases})*100"
    )
    logging.info(
        f"Percentage of accurately called bases is {accurate_base_calling_percentage*100}%"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
