import json
import os
from datetime import datetime

from tqdm import tqdm

from csmed.csmed_cochrane.prepare_dataset import prepare_dataset

if __name__ == "__main__":
    REPOSITORY_PATH = "../../Systematic_Reviews_Update/"
    OUTPUT_DATA_PATH = "../data/external/sr_updates/"

    if not os.path.exists(OUTPUT_DATA_PATH):
        os.makedirs(OUTPUT_DATA_PATH)

    dataset_split = "test"

    sr_updates_reviews = [
        x.split(".")[0]
        for x in os.listdir(f"{REPOSITORY_PATH}/Reviews-Information/")
        if x.endswith(".keywords.txt")
    ]

    added_reviews = []
    index_id = 1
    unsuccessful_reviews = []

    for cochrane_id in tqdm(sr_updates_reviews):
        sr_updates_dataset = prepare_dataset(
            review_id=cochrane_id, output_data_path=OUTPUT_DATA_PATH
        )
        if sr_updates_dataset == {}:
            print(f"No data found for {cochrane_id}")
            unsuccessful_reviews.append(cochrane_id)
            continue

        sr_updates_dataset["cochrane_id"] = cochrane_id
        sr_updates_dataset["dataset_id"] = cochrane_id
        sr_updates_dataset["used_in"] = ["sr_updates"]
        sr_updates_dataset["dataset_split"] = dataset_split

        added_reviews.append(sr_updates_dataset)
        index_id += 1

        final_data = {
            "data_format_version": "1.0.0",
            "updated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "data": added_reviews,
            "unsuccessful_reviews": unsuccessful_reviews,
        }
        with open(f"{OUTPUT_DATA_PATH}/data_index.json", "w") as f:
            json.dump(final_data, f, indent=2)

    print("Done!")
