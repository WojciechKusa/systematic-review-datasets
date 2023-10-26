import json
import os
from datetime import datetime

from tqdm import tqdm

from ec2s.csmed_cochrane.prepare_dataset import prepare_dataset

if __name__ == "__main__":
    EC2S_REVIEWS_FILE = "../data/raw/ec2s_reviews.json"
    OUTPUT_DATA_PATH = "../data/external/ec2s/"

    if not os.path.exists(OUTPUT_DATA_PATH):
        os.makedirs(OUTPUT_DATA_PATH)

    with open(EC2S_REVIEWS_FILE, "r") as f:
        ec2s_reviews = json.load(f)

    dataset_split = "train"
    added_reviews = []
    index_id = 1

    for review in tqdm(ec2s_reviews["data"]):
        review_id = review["cochrane_id"]

        ec2s_dataset = prepare_dataset(
            review_id=review_id, output_data_path=OUTPUT_DATA_PATH
        )
        ec2s_dataset["id"] = index_id
        ec2s_dataset["used_in"] = ["EC2S"]
        ec2s_dataset["dataset_split"] = dataset_split

        added_reviews.append(ec2s_dataset)
        index_id += 1

        final_data = {
            "data_format_version": "1.0.0",
            "updated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "data": added_reviews,
        }
        with open(f"{OUTPUT_DATA_PATH}/data_index.json", "w") as f:
            json.dump(final_data, f, indent=2)
