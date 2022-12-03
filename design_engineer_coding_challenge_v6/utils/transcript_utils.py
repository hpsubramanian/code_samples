import logging
import sys
from collections import Counter
from itertools import chain

import requests


def get_complete_protein_coding_transcripts(gene_object):
    """Function to get complete protein coding transcripts for a given gene
    Complete and protein coding are two of many annotations available for transcripts.
    Complete refers to the presence of a start and stop codon

    Parameters
    ----------
    gene_object : object
        pyensembl gene object

    Returns
    -------
    complete_protein_coding_transcripts : list
        list of protein coding transcripts for a gene
    """
    complete_protein_coding_transcripts = []
    for transcripts in gene_object.transcripts:
        if transcripts.biotype == "protein_coding" and transcripts.complete:
            complete_protein_coding_transcripts.append(transcripts.id)
    logging.info(
        f"{len(complete_protein_coding_transcripts)} complete protein coding transcripts available for {gene_object.gene_name},"
    )
    return complete_protein_coding_transcripts


def get_exons_for_transcripts(pyensembl_handler, transcript_id_list):
    """Function to get ensembl exon ids in a transcript given a ensembl transcript id

    Parameters
    ----------
    pyensembl_handler : object
        pyensembl data handler to query pyensembl database
    transcript_id_list : list
        list of transcripts for which exon ids are returned

    Returns
    -------
    exon_ids : list
        list of exons in all transcripts in the transcript_id_list
    """
    exon_ids = []
    for transcripts in transcript_id_list:
        exon_ids.append(pyensembl_handler.exon_ids_of_transcript_id(transcripts))
    all_exon_ids = list(chain(*exon_ids))
    return all_exon_ids


def count_exon_usage(exon_id_list):
    """Function to count the number of times an exon is used in ensembl transcripts of a gene

    Parameters
    ----------
    exon_id_list : list
        List of exons to be counted.

    Returns
    -------
    exon_counter : Counter
        A counter object containing exons and the number of times it was seen in the list
    """
    exon_counter = Counter(exon_id_list)
    return exon_counter


def find_top_n_exons(exon_counter, percentage_value=0.2):
    """Function to find the top n% of most used exons for given ensembl gene
    Parameters
    ----------
    exon_counter : Counter
        Counter object contain exons and exon counts
    percentage_value : float
        top % of exons to return

    Returns
    -------
    top_exons : list
        top n% of most used exons for given ensembl gene
    """
    # TODO rename percentage value, value varies form 0-1, decimals
    # TODO make percentage value a input variable
    no_of_exons_to_return = round(len(exon_counter) * percentage_value)
    top_exons = [
        exon_id for exon_id, count in exon_counter.most_common(no_of_exons_to_return)
    ]
    return top_exons


def get_canonical_transcript(ensembl_gene_id):
    """Function to get canonical transcript of a given ensembl gene

    Parameters
    ----------
    ensembl_gene_id : str
        ensembl gene id
    Returns
    -------
    canonical_transcript : str
        canonical_transcript of a ensembl gene
    """
    # See https://rest.ensembl.org/documentation/info/lookup
    server = "http://rest.ensembl.org"  # TODO Specify this in a config file.
    ext = f"/lookup/id/{ensembl_gene_id}?expand=0"

    response = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not response.ok:
        response.raise_for_status()
        sys.exit()

    decoded_response = response.json()
    canonical_transcript, version = decoded_response["canonical_transcript"].split(".")
    return canonical_transcript
