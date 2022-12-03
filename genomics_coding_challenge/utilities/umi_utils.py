import logging
from collections import Counter, defaultdict
import re

import pandas as pd
import numpy as np
from scipy import spatial
from Levenshtein import hamming
import pyfastx


def extract_umi_from_read(read, umi_prefix, umi_suffix):
    """Extract UMI from each read using regex

    Args:
        read (str): Read sequence
        umi_prefix (str): UMI prefix
        umi_suffix (str): UMI suffix

    Returns:
        str: UMI sequence
    """
    umi_regex_pattern = re.compile(f"{umi_prefix}(.*){umi_suffix}", re.IGNORECASE)
    umi = umi_regex_pattern.search(read).group(1)
    return umi


def extract_umi_from_fastq(input_fastq, umi_prefix, umi_suffix):
    """Extract all UMI's between a umi prefix and suffix, given a fastq file. Can be gziped.

    Args:
        input_fastq (file): *.fastq or *.fastq.gz file.
        umi_prefix (str): UMI prefix
        umi_suffix (str):  UMI suffix

    Returns:
        pandas dataframe: main_readname_umi_df - A dataframe that contains all sequence headers and UMI sequences in the corresponding sequence
    """
    logging.info("Parsing FASTQ file and extracting UMI's from reads")
    main_readname_umi_df = pd.DataFrame(columns=["ReadName", "UMI"])
    fastq_iterator = pyfastx.Fastx(input_fastq)

    # Ignore quality and comment for now
    for name, sequence, *_ in fastq_iterator:
        umi = extract_umi_from_read(sequence, umi_prefix, umi_suffix)
        current_readname_umi_df = pd.DataFrame({"ReadName": [name], "UMI": [umi]})
        main_readname_umi_df = pd.concat(
            [main_readname_umi_df, current_readname_umi_df], ignore_index=True
        )

    return main_readname_umi_df


def write_data_frame_to_file(dataframe, filename, write_index=None):
    """Function to save a dataframe to a file.

    Args:
        dataframe (pandas df): Pandas dataframe to save to file
        filename (str): Path to write the dataframe to
        write_index (bool, optional): Write row headers if True. Defaults to None.
    """
    logging.info(f"Writing dataframe to file {filename}")
    dataframe.to_csv(filename, header=True, index=write_index, sep="\t", mode="w")


def calculate_pairwise_hamming_distance(unique_umi_counts):
    """Function to calculate pairwise hamming distance between UMI sequences

    Args:
        unique_umi_counts (Counter): Counter object that contains UMI sequences and counts

    Returns:
        pdist_distance_matrix: condensed hamming distance matrix
        squared_distance_matrix: full hamming distance matrix (values across both sides of the diagnol)
    """
    logging.info("Calculating pairwise hamming distance")

    # Convert to array
    transformed_umi_array = np.array(sorted(unique_umi_counts.keys())).reshape(-1, 1)

    # Calculate Distance matrix (condensed) using pdist from scipy
    # Hamming takes care of sequence length
    pdist_distance_matrix = spatial.distance.pdist(
        transformed_umi_array, lambda x, y: hamming(x[0], y[0])
    )

    # Convert it to full matrix
    squared_distance_matrix = pd.DataFrame(
        spatial.distance.squareform(pdist_distance_matrix)
    )

    # Apply row names and column names
    squared_distance_matrix.columns = pd.unique(sorted(unique_umi_counts.keys()))
    squared_distance_matrix.index = pd.unique(sorted(unique_umi_counts.keys()))

    return pdist_distance_matrix, squared_distance_matrix


def get_unique_umi_counts(umi_dataframe):
    """Function applies Counter function to list of all UMIs found in the dataframe.

    Args:
        umi_dataframe (df): Pandas dataframe containing sequence header and UMI

    Returns:
        Counter: Counter dict containing UMI sequences and their counts in the FASTQ file
    """
    logging.info("Obtaining unique UMI count")
    return Counter(umi_dataframe.UMI.tolist())


def calculate_total_umi_bases_sequenced(unique_umi_counts):
    """Function to calculate total UMI bases sequenced

    Args:
        unique_umi_counts (Counter): Counter dict containing UMI sequences and their counts in the FASTQ file

    Returns:
        total_umi_bases: Total number of UMI bases in the fastq file
    """
    total_umi_bases = 0
    for sequences in sorted(unique_umi_counts.keys()):
        total_umi_bases += len(sequences) * unique_umi_counts[sequences]
    return total_umi_bases


def find_similar_umis(hamming_distance_matrix, distance_value):
    """Function to find UMI sequences within a hamming distance value specified.

    Args:
        hamming_distance_matrix (pandas dataframe): Full pandas dataframe containing hamming distance values
        distance_value (int): Find sequences within this hamming distance value

    Returns:
        deduplicated_similar_umis: a list of lists containing similar umi's within the hamming distance threshold
    """
    all_similar_umis = defaultdict(list)

    for columns in hamming_distance_matrix.columns:
        for rows in hamming_distance_matrix.index:
            hamming_distance = hamming_distance_matrix[columns][rows]
            if hamming_distance != 0 and hamming_distance <= distance_value:

                # Explore alternatives, avoid unnecssary looping everytime.
                if columns not in all_similar_umis[columns]:
                    all_similar_umis[columns].append(columns)

                # Explore alternatives, avoid unnecssary looping everytime.
                if rows not in all_similar_umis[columns]:
                    all_similar_umis[columns].append(rows)

                # Explore alternatives, avoid unnecssary sorting everytime.
                all_similar_umis[columns] = sorted(all_similar_umis[columns])

    deduplicated_similar_umis = []

    for similar_umis in all_similar_umis.values():
        if sorted(similar_umis) not in deduplicated_similar_umis:
            deduplicated_similar_umis.append(sorted(similar_umis))

    return deduplicated_similar_umis


def collapse_umis_find_error_bases(unique_umi_counts, similar_umis_list_of_lists):
    """Function to calculate UMI error bases in the fastq file, using information from find_similar_umis

    Args:
        unique_umi_counts (Counter): Counter object that contains UMI sequences and counts
        similar_umis_list_of_lists (list): a list of lists containing similar umi's within the hamming distance threshold

    Returns:
        total_number_sequencing_errors: Total number of UMI bases that were incorrectly called.
    """

    total_number_sequencing_errors = 0

    for umi_list in similar_umis_list_of_lists:

        correct_index = ""
        correct_index_count = 0

        for each_umi in umi_list:
            if unique_umi_counts[each_umi] > correct_index_count:
                correct_index_count = unique_umi_counts[each_umi]
                correct_index = each_umi

        for each_umi in umi_list:
            distance_hamming = hamming(correct_index, each_umi)
            total_number_sequencing_errors += (
                distance_hamming * unique_umi_counts[each_umi]
            )

    return total_number_sequencing_errors
