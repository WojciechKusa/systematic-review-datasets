"""This module prepares TAR 2019 systematic review data."""
import copy
import http.client
import json
import re

import pandas as pd
from Bio import Entrez, Medline

from ec2s.big_screening.utils import (
    set_entrez_email,
    is_prepared,
    save_checksum,
    mark_all_files_prepared,
)

set_entrez_email()

cochrane_id_pattern = r"CD(\d+)"


def get_from_pubmed(df) -> pd.DataFrame:
    """
    :param df: input dataframe containing (PubMedId, Label) dataset.
    """
    labels_column: str = "Label"
    pubmed_id_column: str = "PMID"

    data_size = len(df)
    docs = []
    step = 2000
    for x in range(0, data_size, step):
        pubmed_id_list = df[pubmed_id_column].tolist()[x : x + step]
        handle = Entrez.efetch(
            db="pubmed", id=pubmed_id_list, rettype="medline", retmode="text"
        )

        articles = Medline.parse(handle)
        try:
            articles = list(articles)
        except http.client.HTTPException as e:
            print(e)
            print(pubmed_id_list)
            continue

        for article in articles:
            docs.append(article)

    output_df = pd.DataFrame(docs)
    output_df = output_df.rename(columns={"TI": "Title", "AB": "Abstract"})

    output_df[labels_column] = df[
        df[pubmed_id_column].isin(output_df["PMID"].astype(int).tolist())
    ][labels_column].tolist()

    return output_df


def temporal_search_pubmed(
    query: str, start_date: str, end_date: str
) -> tuple[int, list[str]]:
    if query == "":
        return 0, []
    if start_date:
        start_date = start_date.replace("/", "-")
    if end_date:
        end_date = end_date.replace("/", "-")

    handle = Entrez.esearch(
        db="pubmed", term=query, retmax=1000, mindate=start_date, maxdate=end_date
    )
    try:
        record = Entrez.read(handle)
    except (http.client.HTTPException, RuntimeError) as e:
        print(e)
        return 0, []
    count = int(record["Count"])
    id_list = record["IdList"]
    return count, id_list


def prepare_dataset(input_folder: str, output_folder: str) -> None:
    if is_prepared(output_folder):
        print("PubMed data is already prepared.")
        return

    qrels_df = pd.read_csv(
        f"{input_folder}/sigir2017.qrels",
        sep="\s+",
        header=None,
        names=["review_id", "0", "PMID", "Label"],
    )
    with open(f"{input_folder}/queries_unannotated.json") as f:
        queries = json.load(f)
    with open(f"{input_folder}/systematic_reviews.json") as f:
        reviews_mapping = json.load(f)

    print("PubMed data is being downloaded. This may take a while for the first time.")
    for review_id in qrels_df["review_id"].unique():  # review_ids are just numbers

        review_url = [r["url"] for r in reviews_mapping if r["id"] == review_id][0]
        match = re.search(cochrane_id_pattern, review_url)
        if match:
            cochrane_id = match.group(0)
        else:
            raise ValueError("Review ID not found")

        query = [q for q in queries if q["document_id"] == review_id]
        if len(query) > 0:
            n_found, id_list = temporal_search_pubmed(
                query[0]["query"], query[0]["start_date"], query[0]["end_date"]
            )
        else:
            n_found = 0
            id_list = []

        review_df = copy.copy(qrels_df[qrels_df["review_id"] == review_id])
        print(f"{cochrane_id=}, {len(review_df)=}, {len(id_list)=}")
        if len(id_list) > 0:
            review_df = pd.concat(
                [
                    review_df,
                    pd.DataFrame(
                        {
                            "review_id": review_id,
                            "0": 0,
                            "PMID": id_list,
                            "Label": 0,
                        }
                    ),
                ],
                ignore_index=True,
            )
        review_df["PMID"] = review_df["PMID"].astype(int)
        review_df = review_df.drop_duplicates(subset=["PMID"])

        review_df = get_from_pubmed(review_df)

        review_df["review_id"] = cochrane_id
        review_df["Label"] = review_df["Label"].apply(
            lambda x: 0 if x in [0, 1, 3] else 1
        )

        review_df.to_csv(f"{output_folder}/{cochrane_id}.csv", index=False)
        print(f"Prepared review size: {len(review_df)}")
        save_checksum(
            file=f"{output_folder}/{cochrane_id}.csv",
            dataset_directory=output_folder,
        )

    mark_all_files_prepared(output_folder)
