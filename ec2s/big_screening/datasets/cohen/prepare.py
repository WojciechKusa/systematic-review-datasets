"""This module prepares Drug systematic review data.
Inout data comes from: https://dmice.ohsu.edu/cohenaa/epc-ir-data/epc-ir.clean.tsv
"""
import argparse
import copy

import pandas as pd
from Bio import Entrez, Medline

from ec2s.big_screening.utils import (
    set_entrez_email,
    is_prepared,
    save_checksum,
    mark_all_files_prepared,
)

set_entrez_email()

DATA_URL = "https://dmice.ohsu.edu/cohenaa/epc-ir-data/epc-ir.clean.tsv"

REVIEWS = (
    pd.read_csv(
        DATA_URL,
        sep="\t",
        names=[
            "Drug review topic",
            "EndNote ID",
            "PubMed ID",
            "Abstract Triage Status",
            "Article Triage Status",
        ],
    )["Drug review topic"]
    .unique()
    .tolist()
)


def prepare_dataset(
        output_folder: str,
) -> None:
    if is_prepared(output_folder):
        return

    labels_column: str = "Label"
    pubmed_id_column: str = "PubMed ID"

    df = pd.read_csv(
        DATA_URL,
        sep="\t",
        names=[
            "Drug review topic",
            "EndNote ID",
            "PubMed ID",
            "Abstract Triage Status",
            "Article Triage Status",
        ],
    )

    for review in df["Drug review topic"].unique().tolist():
        review_df = copy.copy(df[df["Drug review topic"] == review])
        print(f"{review=}, {len(review_df)=}")

        review_df[labels_column] = 0
        review_df.loc[review_df["Article Triage Status"] == "I", labels_column] = 1

        pubmed_id_list = review_df[pubmed_id_column].tolist()
        handle = Entrez.efetch(
            db="pubmed", id=pubmed_id_list, rettype="medline", retmode="text"
        )
        articles = Medline.parse(handle)

        docs = []
        for article in articles:
            docs.append(article)
        output_df = pd.DataFrame(docs)
        output_df = output_df.rename(columns={"TI": "Title", "AB": "Abstract"})

        output_df[labels_column] = review_df[
            review_df["PubMed ID"].isin(output_df["PMID"].astype(int).tolist())
        ][labels_column].tolist()

        outfile = f"{output_folder}/{review.strip()}.tsv"
        output_df.to_csv(outfile, sep="\t", index=False)
        save_checksum(file=outfile, dataset_directory=output_folder)

    mark_all_files_prepared(dataset_directory=output_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_folder",
        default="../../data/processed/",
        type=str,
    )

    args = parser.parse_args()

    prepare_dataset(output_folder=args.output_folder)
