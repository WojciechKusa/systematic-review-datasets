import hashlib
import json
import logging
import os
from typing import Any, Optional

import pandas as pd
import requests
from bs4 import BeautifulSoup

from csmed.csmed_cochrane.cochrane_web_parser import (
    parse_cochrane_references,
    parse_data_and_analyses_section,
    parse_search_strategy,
    parse_eligibility_criteria,
    get_safe_soup,
)
from csmed.csmed_cochrane.match_references import expand_references_details

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# get the cookie from your browser and paste it into this file
with open("../data/cookie.txt", "r", encoding="utf-8") as f:
    cookie = f.read().strip()

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/108.0.0.0 Safari/537.36",
    "Cookie": cookie,
}


def _get_versions(review_id: str) -> dict[int, str]:
    _url = f"https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.{review_id}/information"
    soup: Optional[BeautifulSoup] = get_safe_soup(url=_url, headers=HEADERS)
    if not soup:
        return {}

    version_history = soup.find("section", {"class": "versionHistory"})
    versions = [
        version_row.findAll("td")[4].text.strip()
        for version_row in version_history.find("tbody").findAll("tr", recursive=False)
        if version_row.find("td").attrs.get("class") != ["minor-version-table"]
    ]
    return {index_i + 1: version for index_i, version in enumerate(reversed(versions))}


def _get_file(url: str, headers: dict[str, str], output_file: str) -> str:
    """
    Download a file from URL and save it to a file. Return the md5 hash of the file.
    :param url:
    :param headers:
    :param output_file:
    :return: hash of the file, or empty string if the file could not be downloaded
    """
    try:
        r = requests.get(url, headers=headers)
    except requests.exceptions.TooManyRedirects:
        logger.error(f"Too many redirects for {url}")
        return ""
    if r.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(r.content)
        logger.debug(f"Downloaded {url} to {output_file}")
        return hashlib.md5(r.content).hexdigest()
    else:
        logger.warning(f"{r.status_code} while downloading {url} file")
        return ""


def _get_most_recent_available_version(
    review_id: str, review_versions: dict[int, str]
) -> int:
    """
    Check which review version is the most recent one available on Cochrane Library.
    :param review_id:
    :param review_versions:
    :return:
    """
    for version, version_name in reversed(review_versions.items()):
        cochrane_home = f"https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.{review_id}.pub{version}"
        cochrane_revman = f"{cochrane_home}/media/CDSR/{review_id}/table_n/{review_id}StatsDataOnly.rm5"
        try:
            r = requests.get(cochrane_revman, headers=HEADERS)
            if r.status_code == 200 and not r.content.decode().startswith(
                "<!DOCTYPE html>"
            ):
                return version
        except requests.exceptions.TooManyRedirects:
            logger.warning(f"Too many redirects for {cochrane_revman}")
            return version


def _get_cochrane_home_soup(review_url: str) -> BeautifulSoup:
    """Get the soup of the Cochrane Library page of the review.

    :param review_url:
    :return: BeautifulSoup object
    """
    r = requests.get(review_url, headers=HEADERS)
    soup = BeautifulSoup(r.text, "html.parser")
    return soup


def _get_review_title(soup: BeautifulSoup) -> str:
    return soup.find("h1", {"class": "publication-title"}).text.strip()


def _get_review_type(soup: BeautifulSoup) -> str:
    review_type = soup.find("span", {"class": "publish-type"}).text.strip()
    return review_type.split(" - ")[-1]


def _get_review_abstract(soup: BeautifulSoup) -> str:
    """Returns the abstract of the review from Cochrane Library as a raw text."""
    return (
        soup.find("section", {"class": "abstract"})
        .find("div", {"class": "full_abstract"})
        .text.strip()
    )


def get_review_doi(review_id: str, review_version: int) -> str:
    if review_version == 1:
        return f"https://doi.org/10.1002/14651858.{review_id}"
    else:
        return f"https://doi.org/10.1002/14651858.{review_id}.pub{review_version}"


def prepare_dataset(
    review_id: str, output_data_path: str, use_version: str = "latest"
) -> dict[str, Any]:
    """This script downloads the pdf and revman from Cochrane Library and saves it to output_data_path.
    It also parses the data and returns a dict with data.

    :param use_version: "latest" or "all" or numerical version as a string, e.g. "1"
    :param review_id: Cochrane review ID, e.g. "CD010771"
    :param output_data_path: path to the output directory where the data will be saved
    :return:
    """
    if use_version not in ["latest", "all"] and not use_version.isdigit():
        raise ValueError(
            f"Wrong value for use_version: {use_version}. Should be latest, all or a version number."
        )

    versions = _get_versions(review_id=review_id)
    if use_version == "latest":
        if len(versions) == 0:
            logger.error(f"Review {review_id} has no versions")
            return {}
        if list(versions.keys())[-1] != 1:
            logger.debug(
                f"Review {review_id} has more than one version. Checking latest available version"
            )
            version = _get_most_recent_available_version(review_id, versions)
            if not version:
                logger.error(f"Review {review_id} has no available versions with RevMan files")
                return {}

            logger.debug(
                f"Most recent available version is {version} out of {list(versions.keys())[-1]}"
            )
            #  there were problems with cochrane website on Feb 23
            if review_id == "CD010038":
                version = 2
        else:
            version = 1
    elif use_version == "all":
        if len(versions) == 0:
            logger.warning(f"Review {review_id} has no versions")
            return {}
        version = "all"
        raise NotImplementedError("Downloading all versions is not implemented yet")
    else:
        version = int(use_version)

    if version == 1:
        logger.debug(f"Review {review_id} has only one version")
        cochrane_home = (
            f"https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.{review_id}"
        )
        cochrane_pdf = (
            f"{cochrane_home}/pdf/CDSR/{review_id}/rel0001/{review_id}/{review_id}.pdf"
        )
    else:
        cochrane_home = f"https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.{review_id}.pub{version}"
        cochrane_pdf = f"{cochrane_home}/pdf/CDSR/{review_id}/{review_id}.pdf"

    doi = get_review_doi(review_id=review_id, review_version=version)
    cochrane_revman = (
        f"{cochrane_home}/media/CDSR/{review_id}/table_n/{review_id}StatsDataOnly.rm5"
    )

    out_path = f"{output_data_path}/{review_id}/"
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    pdf_md5 = _get_file(
        url=cochrane_pdf, headers=HEADERS, output_file=f"{out_path}/{review_id}.pdf"
    )
    revman_md5 = _get_file(
        url=cochrane_revman, headers=HEADERS, output_file=f"{out_path}/{review_id}.rm5"
    )

    cochrane_home_soup: Optional[BeautifulSoup] = get_safe_soup(
        url=cochrane_home, headers=HEADERS
    )
    if not cochrane_home_soup:
        logger.error(
            f"Review {review_id} has no Cochrane home page or it is not available"
        )
        return {}

    title = _get_review_title(soup=cochrane_home_soup)
    abstract = _get_review_abstract(soup=cochrane_home_soup)
    review_type = _get_review_type(soup=cochrane_home_soup)
    criteria = parse_eligibility_criteria(soup=cochrane_home_soup)

    cochrane_appendix = f"{cochrane_home}/appendices"
    cochrane_appendix_soup: Optional[BeautifulSoup] = get_safe_soup(
        url=cochrane_appendix, headers=HEADERS
    )
    if cochrane_appendix_soup:
        search_strategy = parse_search_strategy(soup=cochrane_appendix_soup)
    else:
        search_strategy = ""

    cochrane_references = f"{cochrane_home}/references"
    cochrane_references_soup: Optional[BeautifulSoup] = get_safe_soup(
        url=cochrane_references, headers=HEADERS
    )
    df = parse_cochrane_references(soup=cochrane_references_soup)
    df = expand_references_details(df)
    df.to_csv(f"{out_path}/references.csv", index=False)

    n_studies_included = len(
        df[df["reference_type"] == "included"]["study_id"].unique()
    )
    n_references_included = len(df[df["reference_type"] == "included"])
    unique_references_included = len(
        df[df["reference_type"] == "included"]["citation"].unique()
    )
    unique_references_excluded = len(
        df[df["reference_type"] == "excluded"]["citation"].unique()
    )

    data_df = parse_data_and_analyses_section(soup=cochrane_references_soup, review_id=review_id)
    data_df.to_csv(f"{out_path}/data_and_analyses.csv", index=False)

    if not data_df.empty:
        n_comparisons = len(data_df["comparison_name"].unique())
        try:
            n_outcomes = len(data_df["outcome_name"].unique())
            n_outcomes_and_subgroups = len(data_df)
            avg_studies_for_outcome = (
                data_df[
                    pd.to_numeric(data_df["No. of studies"], errors="coerce").notnull()
                ]["No. of studies"]
                .astype(int)
                .mean()
            )
        except (KeyError, ValueError):  # fixme: check other option for DTA reviews
            n_outcomes = 0
            n_outcomes_and_subgroups = 0
            avg_studies_for_outcome = 0
    else:
        n_comparisons = 0
        n_outcomes = 0
        n_outcomes_and_subgroups = 0
        avg_studies_for_outcome = 0

    review_details = {
        "title": title,
        "abstract": abstract,
        "review_type": review_type,
        "doi": doi,
        "review_id": review_id,
        "criteria": criteria,
        "search_strategy": search_strategy,
    }
    with open(f"{out_path}/review_details.json", "w") as f:
        json.dump(review_details, f, indent=4)

    return {
        "title": title,
        "doi": doi,
        "review_id": review_id,
        "source": "Cochrane",
        "versions": versions,
        "url": cochrane_home,
        "pdf": cochrane_pdf,
        "pdf_md5": pdf_md5,
        "revman": cochrane_revman,
        "revman_md5": revman_md5,
        "references": cochrane_references,
        "n_studies_included": n_studies_included,
        "n_references_included": n_references_included,
        "unique_references_included": unique_references_included,
        "unique_references_excluded": unique_references_excluded,
        "n_comparisons": n_comparisons,
        "n_outcomes": n_outcomes,
        "n_outcomes_and_subgroups": n_outcomes_and_subgroups,
        "avg_studies_for_outcome": avg_studies_for_outcome,
    }
