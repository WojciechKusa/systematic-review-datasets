import json
import os
from datetime import datetime

from tqdm import tqdm

from ec2s.ec2s.prepare_dataset import prepare_dataset

if __name__ == "__main__":
    REVIEWS_FILE = "../data/raw/pcs_reviews.json"
    OUTPUT_DATA_PATH = "../data/external/pcs/"

    if not os.path.exists(OUTPUT_DATA_PATH):
        os.makedirs(OUTPUT_DATA_PATH)

    dataset_split = "Testing"

    with open(REVIEWS_FILE, "r") as f:
        prospective_reviews = json.load(f)

    added_reviews = []

    for index_id, review in tqdm(enumerate(prospective_reviews["data"])):
        review_id = review["cochrane_id"]

        pcs_dataset = prepare_dataset(
            review_id=review_id, output_data_path=OUTPUT_DATA_PATH, use_version="2"
        )
        pcs_dataset["cochrane_id"] = review_id
        pcs_dataset["dataset_id"] = index_id
        pcs_dataset["used_in"] = ["pcs"]
        pcs_dataset["dataset_split"] = dataset_split

        added_reviews.append(pcs_dataset)

        final_data = {
            "data_format_version": "1.0.0",
            "updated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "data": added_reviews,
        }
        with open(f"{OUTPUT_DATA_PATH}/data_index.json", "w") as f:
            json.dump(final_data, f, indent=2)
