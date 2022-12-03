import logging

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import cluster


def create_umi_frequency_distribution_graph(unique_umi_counts, output_file_name):
    """This function plots UMI count Frequency distribution. (X-axis: Number of times a barcode was observed in the library; Y-axis: Number of barcodes observed that many times)
    It also returns the data used to plot as a table.

    Args:
        unique_umi_counts (Counter Object): Counter object containing UMI sequence and corresponding counts
        output_file_name (Plot Image): File path to plot UMI Frequency distribution.

    Returns:
        freq_df: Pandas dataframe, containing UMI count frequency distribution
    """

    logging.info("Plotting UMI Count Frequency distribution")

    umi_counts_df = pd.DataFrame({"Frequency": sorted(unique_umi_counts.values())})

    try:
        # Plot and Save
        umi_counts_df["Frequency"].value_counts().plot(
            kind="bar",
            ylabel="No. of UMI's with the same count",
            xlabel="UMI Count in the Library",
        )
        plt.yticks(
            np.arange(
                min(umi_counts_df["Frequency"].value_counts()),
                max(umi_counts_df["Frequency"].value_counts()) + 1,
                1.0,
            )
        )
        plt.xticks(rotation=30, ha="right")
        plt.savefig(output_file_name, bbox_inches="tight", dpi=600)

        logging.info(f"UMI Frequency Distribution plot saved to {output_file_name}")
    except Exception:
        raise "Failed to plot UMI Count Frequency"

    # Get data as a pandas df
    freq_df = pd.DataFrame(umi_counts_df["Frequency"].value_counts())

    return freq_df.rename_axis("UMI Counts in Library").reset_index()


def create_dendrogram(distance_matrix, clustering_method, output_file_name, labels):
    """This function plots dendrogram, given a pdist distance matrix (from pairwise hamming distance calculation).

    Args:
        distance_matrix (pdist distance matrix): condensed distance matrix from pairwise hamming distance calculation
        clustering_method (_type_): Clustering method obtained from command line. \
            See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
        output_file_name (str): File path to plot dendrogram
        labels (list): labels for dendrogram
    """

    logging.info("Plotting dendrogram")

    # Performs hierarchical/agglomerative clustering
    linkage = cluster.hierarchy.linkage(distance_matrix, method=clustering_method)

    try:
        # Plot and save dendrogram
        plt.figure()
        plt.title("Hierarchical Clustering Dendrogram")
        plt.xlabel(f"Distance (clustering method - {clustering_method})")
        cluster.hierarchy.dendrogram(
            linkage, orientation="left", labels=labels, show_leaf_counts=True
        )
    except Exception:
        raise "Failed to plot dendrogram"
    plt.savefig(output_file_name, bbox_inches="tight", dpi=600)
    return
