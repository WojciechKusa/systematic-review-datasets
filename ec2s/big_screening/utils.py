import hashlib
import json
import os

from Bio import Entrez


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
    else:
        with open(f"{dataset_directory}/{checksum_file}", "r", encoding="utf-8") as f:
            checksums = json.load(f)

        if ALL_FILES_PREPARED_FLAG not in checksums or not checksums[ALL_FILES_PREPARED_FLAG]:
            return False

        for filename, checksum in checksums.items():
            if filename == ALL_FILES_PREPARED_FLAG:
                continue
            with open(f"{dataset_directory}/{filename}", "rb") as f:
                md5sum = hashlib.md5(f.read()).hexdigest()

            if checksum == md5sum:
                continue
            else:
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


def mark_all_files_prepared(dataset_directory: str, checksum_file: str = "_checksums.json") -> None:
    """Adds additional flag to the checksum file to mark all files as prepared."""
    if not os.path.exists(f"{dataset_directory}/{checksum_file}"):
        checksums = {}
    else:
        with open(f"{dataset_directory}/{checksum_file}", "r", encoding="utf-8") as f:
            checksums = json.load(f)

    checksums[ALL_FILES_PREPARED_FLAG] = True

    with open(f"{dataset_directory}/{checksum_file}", "w", encoding="utf-8") as f:
        json.dump(checksums, f, indent=4)