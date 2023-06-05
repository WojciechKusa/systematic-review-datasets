"""This module prepares TAR 2019 systematic review data."""
import http.client

import pandas as pd
from Bio import Entrez, Medline

from ec2s.big_screening.utils import (
    set_entrez_email,
    is_prepared,
    save_checksum,
    mark_all_files_prepared,
)

set_entrez_email()


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


def prepare_dataset(
    input_folder: str, output_folder: str, dataset_splits: dict[str, dict[str, str]]
) -> None:
    if is_prepared(output_folder):
        return

    for dataset_split, review_types in dataset_splits.items():
        for review_type, qrels_file in review_types.items():
            qrels_df = pd.read_csv(
                f"{input_folder}/tar-master/2019-TAR/Task2/{dataset_split}/{review_type}/qrels/{qrels_file}",
                sep="\s+",
                header=None,
                names=["review_id", "0", "PMID", "Label"],
            )

            print(
                "PubMed data is being downloaded. This may take a while for the first time."
            )
            for review_id in qrels_df["review_id"].unique():
                review_df = qrels_df[qrels_df["review_id"] == review_id]
                print(f"{review_id=}, {len(review_df)=}")
                review_df = get_from_pubmed(review_df)
                review_df.to_csv(f"{output_folder}/{review_id}.csv", index=False)
                print(f"Prepared review size: {len(review_df)}")
                save_checksum(
                    file=f"{output_folder}/{review_id}.csv",
                    dataset_directory=output_folder,
                )

            mark_all_files_prepared(output_folder)
