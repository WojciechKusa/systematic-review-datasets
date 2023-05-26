import logging
import os

import pandas as pd
import requests
from grobid_client.grobid_client import GrobidClient
from tqdm import tqdm

client = GrobidClient(config_path="../data/raw/grobid_config.json")

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def download_file(url, download_path):
    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:55.0) Gecko/20100101 Firefox/55.0",
    }
    r = requests.get(url, stream=True, allow_redirects=True, headers=headers)
    with open(download_path, "wb") as f:
        for ch in r:
            f.write(ch)


if __name__ == "__main__":
    datasets = [
        "../data/external/ec2s/",
        "../data/external/tar-systematic_reviews/",
    ]

    for dataset in datasets:
        reviews = [
            review
            for review in os.listdir(dataset)
            if os.path.isdir(f"{dataset}/{review}")
        ]
        for review in tqdm(reviews):
            print(review)
            df = pd.read_csv(f"{dataset}/{review}/references.csv")
            save_folder = f"{dataset}/{review}/pdfs"

            logger.info("Downloading PDFs")
            for index, row in df.iterrows():
                if "CORE_PDF_LINK" in row:
                    pdf_url = row["CORE_PDF_LINK"]
                elif "S2_PDF_LINK" in row:
                    pdf_url = row["S2_PDF_LINK"]
                else:
                    pdf_url = None

                if not os.path.exists(save_folder):
                    os.makedirs(save_folder)

                if pdf_url and pdf_url == pdf_url:
                    filename = row["reference_id"]
                    download_file(pdf_url, f"{save_folder}/{filename}.pdf")

            logger.info("Processing PDFs")
            client.process("processFulltextDocument", save_folder, n=20, verbose=True)
