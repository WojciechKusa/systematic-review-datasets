import re

import pandas as pd
import requests
from Bio import Entrez, Medline
from tqdm import tqdm
from urllib.error import HTTPError
from http.client import IncompleteRead


S2_API_KEY_FILE = "../data/raw/semanticscholar_api_key.txt"
CORE_API_KEY_FILE = "../data/raw/core_api_key.txt"
S2_API_ENDPOINT = "https://api.semanticscholar.org/graph/v1/paper/"
CORE_API_ENDPOINT = "https://api.core.ac.uk/v3/discover"

try:
    with open(S2_API_KEY_FILE, "r", encoding="utf-8") as f:
        s2_api_key = f.read().strip()
except FileNotFoundError:
    s2_api_key = ""

try:
    with open(CORE_API_KEY_FILE, "r", encoding="utf-8") as f:
        core_api_key = f.read().strip()
except FileNotFoundError:
    core_api_key = ""


def get_paper_id_by_doi(doi):
    url = f"{S2_API_ENDPOINT}DOI:{doi}?fields=url,title,externalIds,citationCount,openAccessPdf"
    r = requests.get(url, headers={"x-api-key": s2_api_key})
    return r.json() if r.status_code == 200 else None


def get_paper_full_text_by_doi_core(doi):
    headers = {
        "Content-Type": "application/json",
        "Authorization": "Bearer " + core_api_key,
    }
    r = requests.post(url=CORE_API_ENDPOINT, headers=headers, json={"doi": doi})
    return r.json() if r.status_code == 200 else None


def search_pubmed(query: str) -> str:
    """Search pubmed for a query and return the pubmed id."""
    if not query:
        return ""
    handle = Entrez.esearch(db="pubmed", term=query)
    try:
        record = Entrez.read(handle)
    except (HTTPError, RuntimeError) as e:
        print(f"Error: {e}")
        return ""
    return str(record["IdList"][0]) if len(record["IdList"]) > 0 else ""


def assign_pubmed_ids(included_df: pd.DataFrame) -> pd.DataFrame:
    """Search for data in pubmed based on the study title and authors."""
    Entrez.email = "test@gmail.com"
    included_df["PUBMED_ID_TITLE_AUTHORS"] = ""
    included_df["PUBMED_ID_TITLE"] = ""
    for index, row in tqdm(included_df.iterrows(), total=len(included_df)):
        title = row["title"]
        authors = row["authors"]

        _pubmed_id = search_pubmed(f"{title} AND {authors}")
        included_df.loc[index, "PUBMED_ID_TITLE_AUTHORS"] = _pubmed_id

        _pubmed_id = search_pubmed(f"{title}")
        if _pubmed_id:
            try:
                handle = Entrez.efetch(
                    db="pubmed", id=_pubmed_id, rettype="medline", retmode="text"
                )
                if record := list(Medline.parse(handle))[0]:
                    _title = str(record.get("TI", "").strip(' .,;:!?-"'))
                    _title = re.sub(r"[^a-zA-Z0-9]+", " ", _title).lower().strip()
                    if (
                        _title
                        != re.sub(r"[^a-zA-Z0-9]+", " ", str(title)).lower().strip()
                    ):
                        _pubmed_id = ""
            except (ValueError, HTTPError, IncompleteRead):
                _pubmed_id = ""
        included_df.loc[index, "PUBMED_ID_TITLE"] = _pubmed_id

        s2 = get_paper_id_by_doi(
            str(row.get("doi", "")).replace("https://doi.org/", "")
        )
        included_df.loc[index, "PUBMED_ID_S2"] = (
            s2["externalIds"].get("PubMed", "") if s2 else ""
        )
        included_df.loc[index, "N_CITATIONS"] = s2["citationCount"] if s2 else ""
        included_df.loc[index, "S2_PDF_LINK"] = (
            s2["openAccessPdf"]["url"] if s2 and s2["openAccessPdf"] else ""
        )

        core = get_paper_full_text_by_doi_core(
            str(row.get("doi", "")).replace("https://doi.org/", "")
        )
        included_df.loc[index, "CORE_PDF_LINK"] = (
            core["fullTextLink"] if core and core["fullTextLink"] else ""
        )

    included_df = normalise_pubmed_ids(df=included_df)
    return included_df


def normalise_pubmed_ids(df: pd.DataFrame) -> pd.DataFrame:
    """There are 4 different columns with pubmed ids. We normalise them to one column.
    Priority: pubmed_id, PUBMED_ID_S2, PUBMED_ID_TITLE_AUTHORS, PUBMED_ID_TITLE
    cochrane web > doi from s2 > title + authors from pubmed > title from pubmed"""
    columns = [
        "pubmed_id",
        "PUBMED_ID_S2",
        "PUBMED_ID_TITLE_AUTHORS",
        "PUBMED_ID_TITLE",
    ]
    for column in columns:
        df[column] = df[column].fillna("")
        df[column] = df[column].astype(str)
        df[column] = df[column].replace("", None)

    df["PUBMED_ID"] = df["pubmed_id"]
    df["PUBMED_ID"] = df["PUBMED_ID"].fillna(df["PUBMED_ID_S2"])
    df["PUBMED_ID"] = df["PUBMED_ID"].fillna(df["PUBMED_ID_TITLE_AUTHORS"])
    df["PUBMED_ID"] = df["PUBMED_ID"].fillna(df["PUBMED_ID_TITLE"])

    df["PUBMED_ID"] = df["PUBMED_ID"].fillna("")
    df["PUBMED_ID"] = df["PUBMED_ID"].astype(str)
    return df


def expand_references_details(web_references_df: pd.DataFrame) -> pd.DataFrame:
    web_references_df = assign_pubmed_ids(web_references_df)
    return web_references_df


def merge_papers_to_studies(
    web_references_df: pd.DataFrame,
    revman_df: pd.DataFrame,
    review_id: str,
    output_path: str,
) -> pd.DataFrame:
    """Merge the papers to the studies."""

    web_references_df = web_references_df[
        web_references_df["reference_type"] == "included"
    ].reset_index()

    merged_df = match_web_to_revman_studies(web_references_df, revman_df)
    not_found_df = web_references_df[
        ~web_references_df["header"].isin(merged_df["header"])
    ]
    not_found_other = revman_df[~revman_df["STUDY_NAME"].isin(merged_df["STUDY_NAME"])]

    if not not_found_df.empty or not not_found_other.empty:
        print(review_id)
        print(not_found_df)
        print(not_found_other)

    print(
        f"{review_id}\t {len(web_references_df)=}, {len(revman_df)=}, {len(merged_df)=}, {len(not_found_df)=}, {len(not_found_other)=}"
    )
    if not not_found_df.empty:
        not_found_df.to_csv(
            f"{output_path}/{review_id}/not_found.csv",
            index=False,
        )

    return assign_pubmed_ids(merged_df)


def match_web_to_revman_studies(
    web_references_df: pd.DataFrame, revman_df: pd.DataFrame
) -> pd.DataFrame:
    """In web_references_df we match column 'header' to the column STUDY_NAME in revman_df.
    First the web_references_df is preprocessed to strip of additional characters."""
    web_references_column = "header"
    revman_column = "STUDY_NAME"

    web_references_df[web_references_column] = web_references_df[
        web_references_column
    ].str.replace(r"\[.*\]", "", regex=True)
    web_references_df[web_references_column] = web_references_df[
        web_references_column
    ].str.replace(r"\{.*\}", "", regex=True)
    web_references_df = normalise_column(web_references_df, web_references_column)

    revman_df = normalise_column(revman_df, revman_column)
    return pd.merge(
        web_references_df,
        revman_df,
        left_on=web_references_column,
        right_on=revman_column,
        how="inner",
    )


def normalise_column(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """Normalise the column by removing brackets and extra spaces.
    Also replace the unicode character '‐' with a hyphen."""
    df[column] = df[column].str.replace(r"\(.*\)", "", regex=True)
    df[column] = df[column].str.replace(r"\s+", " ", regex=True)
    df[column] = df[column].str.replace("‐", "-")
    df[column] = df[column].str.strip()
    return df
