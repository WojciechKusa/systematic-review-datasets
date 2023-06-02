"""This module prepares SWIFT systematic review data."""
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

    df[labels_column] = 1
    df.loc[df["Status"] == "Excluded", labels_column] = 0

    data_size = len(df)
    docs = []
    for x in range(0, data_size, 5000):
        pubmed_id_list = df[pubmed_id_column].tolist()[x : x + 5000]
        handle = Entrez.efetch(
            db="pubmed", id=pubmed_id_list, rettype="medline", retmode="text"
        )

        articles = Medline.parse(handle)
        for article in articles:
            try:
                docs.append(article)
            except http.client.HTTPException as e:
                print(e)
                print(article)

    output_df = pd.DataFrame(docs)
    output_df = output_df.rename(columns={"TI": "Title", "AB": "Abstract"})

    output_df[labels_column] = df[
        df[pubmed_id_column].isin(output_df["PMID"].astype(int).tolist())
    ][labels_column].tolist()

    return output_df


def prepare_fluoride_dataset(df: pd.DataFrame) -> pd.DataFrame:
    """
    :param df: input dataframe containing Fluoride dataset
    """
    labels_column: str = "Label"
    df["Title"] = df["Title"].fillna("")
    df["Abstract"] = df["Abstract"].fillna("")

    df[labels_column] = 1
    df.loc[df["Included"] == "EXCLUDED", labels_column] = 0

    return df


def prepare_neuropathic_pain_dataset(df: pd.DataFrame) -> pd.DataFrame:
    """
    :param df: input dataframe containing NeuropathicPain dataset
    """
    labels_column: str = "Label"
    df["Title"] = df["Title"].fillna("")
    df["Abstract"] = df["Abstract"].fillna("")

    df["tmp_label"] = df["Label"]
    df[labels_column] = 1
    df.loc[df["tmp_label"] == "Excluded", labels_column] = 0

    return df


reviews_version = {
    "Neuropain": prepare_neuropathic_pain_dataset,
    "Fluoride": prepare_fluoride_dataset,
    "BPA": get_from_pubmed,
    "Transgenerational": get_from_pubmed,
    "PFOS-PFOA": get_from_pubmed,
}

file_to_review_mapping = {
    "ohat": [
        "BPA",
        "Transgenerational",
        "PFOS-PFOA",
        "Fluoride",
    ],
    "camrades": ["Neuropain"],
}

REVIEWS = [x.split(".")[0] for x in reviews_version.keys()]


def prepare_dataset(
    input_files: dict[str, str],
    output_folder: str,
) -> None:
    if is_prepared(output_folder):
        return

    print("PubMed data is being downloaded. This may take a while for the first time.")
    for file_name, input_file in input_files.items():
        for review in file_to_review_mapping[file_name]:
            df = pd.read_excel(input_file, sheet_name=review)
            print(f"Processing {review}, {len(df)=}")
            df = reviews_version[review](df)
            df.to_csv(f"{output_folder}/{review}.csv", index=False)
            save_checksum(
                file=f"{output_folder}/{review}.csv", dataset_directory=output_folder
            )

    mark_all_files_prepared(output_folder)
