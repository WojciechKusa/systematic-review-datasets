import logging
import os
import ssl

import pandas as pd
import requests
from grobid_client.grobid_client import GrobidClient
from tqdm import tqdm

client = GrobidClient(config_path="../data/raw/grobid_config.json")

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def download_file(url: str, download_path: str):
    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:55.0) Gecko/20100101 Firefox/55.0",
    }
    try:
        r = requests.get(
            url,
            stream=True,
            allow_redirects=True,
            headers=headers,
            verify=ssl.CERT_NONE,
        )
    except (
        requests.exceptions.InvalidSchema,
        requests.exceptions.ConnectTimeout,
        requests.exceptions.ConnectionError,
        requests.exceptions.ReadTimeout,
        requests.exceptions.TooManyRedirects,
    ) as e:
        logger.warning(f"Invalid schema for {url}, error: {e}")
        return
    with open(download_path, "wb") as f:
        for ch in r:
            f.write(ch)


if __name__ == "__main__":
    datasets = [
        "../data/external/ec2s/",
        "../data/external/pcs/",
    ]

    for dataset in datasets:
        reviews = [
            review
            for review in os.listdir(dataset)
            if os.path.isdir(f"{dataset}/{review}")
        ]
        for review in tqdm(reviews):
            df = pd.read_csv(f"{dataset}/{review}/references.csv")
            save_folder = f"{dataset}/{review}/pdfs"

            for index, row in tqdm(df.iterrows(), total=len(df)):
                if "CORE_PDF_LINK" in row:
                    pdf_url = row["CORE_PDF_LINK"]
                elif "S2_PDF_LINK" in row:
                    pdf_url = row["S2_PDF_LINK"]
                else:
                    pdf_url = None

                if row["reference_type"] not in ["included", "excluded"]:
                    continue

                if not os.path.exists(save_folder):
                    os.makedirs(save_folder)

                if pdf_url and pdf_url == pdf_url:
                    filename = row["reference_id"]
                    download_file(pdf_url, f"{save_folder}/{filename}.pdf")

            client.process("processFulltextDocument", save_folder, n=20, verbose=True)

            extracted_full_texts = []
            for file in os.listdir(save_folder):
                if file.endswith("tei.xml"):
                    try:
                        exclusion_reason = df[df["reference_id"] == file.split(".")[0]][
                            "exclusion_reason"
                        ].values[0]
                    except KeyError:
                        exclusion_reason = ""
                    extracted_full_texts.append(
                        {
                            "review_id": review,
                            "document_id": file.split(".")[0],
                            "reason_for_exclusion": exclusion_reason,
                            "decision": df[df["reference_id"] == file.split(".")[0]][
                                "reference_type"
                            ].values[0],
                        }
                    )

            extracted_full_texts_df = pd.DataFrame(extracted_full_texts)
            if len(extracted_full_texts_df) == 0:
                continue
            extracted_full_texts_df.to_csv(
                f"{dataset}/{review}/extracted_full_texts.csv", index=False
            )
