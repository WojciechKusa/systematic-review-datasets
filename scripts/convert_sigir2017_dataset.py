import itertools
import json
import os
from datetime import datetime
import re

from tqdm import tqdm

from ec2s.ec2s.prepare_dataset import prepare_dataset

review_files = "https://raw.githubusercontent.com/ielab/SIGIR2017-SysRev-Collection/master/systematic_reviews.json"

pattern = r'CD(\d+)'


if __name__ == "__main__":
    REPOSITORY_PATH = "../../SIGIR2017-SysRev-Collection/"
    OUTPUT_DATA_PATH = "../data/external/sigir2017/"

    if not os.path.exists(OUTPUT_DATA_PATH):
        os.makedirs(OUTPUT_DATA_PATH)

    review_type = "Intervention"
    dataset_split = "Training"

    with open(f"{REPOSITORY_PATH}/systematic_reviews.json", "r") as f:
        sigir2017_reviews = json.load(f)

    added_reviews = []
    index_id = 1

    for review in sigir2017_reviews:
        match = re.search(pattern, review["url"])
        if match:
            review_id = match.group(0)
        else:
            raise ValueError("Review ID not found")

        sigir2017_dataset = prepare_dataset(
            review_id=review_id, output_data_path=OUTPUT_DATA_PATH
        )
        sigir2017_dataset["cochrane_id"] = index_id
        sigir2017_dataset["dataset_id"] = review["id"]
        sigir2017_dataset["used_in"] = ["SIGIR2017"]
        sigir2017_dataset["review_type"] = review_type
        sigir2017_dataset["dataset_split"] = "train"

        added_reviews.append(sigir2017_dataset)
        index_id += 1

        final_data = {
            "data_format_version": "1.0.0",
            "updated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "data": added_reviews,
        }
        with open(f"{OUTPUT_DATA_PATH}/data_index.json", "w") as f:
            json.dump(final_data, f, indent=2)
