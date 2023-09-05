import hashlib
import http.client
import json
import os
import time

import pandas as pd
from Bio import Entrez, Medline
from tqdm import tqdm


def set_entrez_email() -> None:
    entrez_email_file = os.path.join(
        os.path.dirname(__file__), "../../config/entrez_email.txt"
    )

    if os.path.exists(entrez_email_file):
        with open(entrez_email_file, "r", encoding="utf-8") as f:
            user_email = f.read().strip()
        Entrez.email = user_email
    elif os.environ.get("ENTREZ_EMAIL"):
        Entrez.email = os.environ.get("ENTREZ_EMAIL")
    else:
        raise ValueError("Please provide an email address for Entrez API")


ALL_FILES_PREPARED_FLAG = "all_files_prepared"


def is_prepared(dataset_directory: str, checksum_file: str = "_checksums.json") -> bool:
    if not os.path.exists(f"{dataset_directory}/{checksum_file}"):
        return False
    with open(f"{dataset_directory}/{checksum_file}", "r", encoding="utf-8") as f:
        checksums = json.load(f)

    if (
        ALL_FILES_PREPARED_FLAG not in checksums
        or not checksums[ALL_FILES_PREPARED_FLAG]
    ):
        return False

    for filename, checksum in checksums.items():
        if filename == ALL_FILES_PREPARED_FLAG:
            continue
        with open(f"{dataset_directory}/{filename}", "rb") as f:
            md5sum = hashlib.md5(f.read()).hexdigest()

        if checksum != md5sum:
            return False

    return True


def save_checksum(
    file: str, dataset_directory: str, checksum_file: str = "_checksums.json"
) -> None:
    """Save the checksum of a file to a JSON file. If the JSON file does not exist, it will be created.
    If the JSON file exists, the checksum will be added to the JSON file.
    If the file already exists in the JSON file, the checksum will be updated."""
    if not os.path.exists(f"{dataset_directory}/{checksum_file}"):
        checksums = {}
    else:
        with open(f"{dataset_directory}/{checksum_file}", "r", encoding="utf-8") as f:
            checksums = json.load(f)

    with open(file, "rb") as f:
        md5sum = hashlib.md5(f.read()).hexdigest()

    checksums[os.path.basename(file)] = md5sum

    with open(f"{dataset_directory}/{checksum_file}", "w", encoding="utf-8") as f:
        json.dump(checksums, f, indent=4)


def mark_all_files_prepared(
    dataset_directory: str, checksum_file: str = "_checksums.json"
) -> None:
    """Adds additional flag to the checksum file to mark all files as prepared."""
    if not os.path.exists(f"{dataset_directory}/{checksum_file}"):
        checksums = {}
    else:
        with open(f"{dataset_directory}/{checksum_file}", "r", encoding="utf-8") as f:
            checksums = json.load(f)

    checksums[ALL_FILES_PREPARED_FLAG] = True

    with open(f"{dataset_directory}/{checksum_file}", "w", encoding="utf-8") as f:
        json.dump(checksums, f, indent=4)


def fetch_pubmed_abstracts(
    pmids: list[str],
    outdir: str,
    batch_size: int = 100,
    delay: float = 0.333,
    overwrite: bool = False,
    verbose: bool = False,
):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not isinstance(pmids, list):
        pmids = list(pmids)

    n_chunks = len(pmids) // batch_size
    for i in tqdm(range(n_chunks + (1 if n_chunks * batch_size < len(pmids) else 0))):
        start, end = i * batch_size, (i + 1) * batch_size
        outfname = f"{outdir}/pubmed.{i}.json"
        if os.path.exists(outfname) and not overwrite:
            continue

        query = ",".join(pmids[start:end])
        handle = Entrez.efetch(
            db="pubmed", id=query, rettype="gb", retmode="xml", retmax=batch_size
        )
        record = Entrez.read(handle)
        if len(record["PubmedArticle"]) != len(pmids[start:end]) and verbose:
            print(
                f"Queried {len(pmids[start:end])}, returned {len(record['PubmedArticle'])}"
            )

        time.sleep(delay)
        # dump to JSON
        with open(outfname, "wt", encoding="utf_8") as file:
            json.dump(record, file, indent=2)


def get_from_pubmed(
    df: pd.DataFrame, labels_column: str = "Label", pubmed_id_column: str = "PMID"
) -> pd.DataFrame:
    """
    :param df: input dataframe containing (PubMedId, Label) dataset.
    :param pubmed_id_column:  column name of PubMedId
    :param labels_column: column name of labels
    """
    set_entrez_email()

    data_size = len(df)
    step = 200
    delay = 20

    docs = []
    for x in tqdm(range(0, data_size, step), unit_scale=step, unit="articles"):
        pubmed_id_list = df[pubmed_id_column].tolist()[x : x + step]
        while True:
            try:
                handle = Entrez.efetch(
                    db="pubmed", id=pubmed_id_list, rettype="medline", retmode="text"
                )

                articles = Medline.parse(handle)
                articles = list(articles)
                break
            except http.client.HTTPException as e:
                print(e)
                print(f"Sleeping for {delay} seconds...")
                time.sleep(delay)

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
    set_entrez_email()

    if not query:
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
