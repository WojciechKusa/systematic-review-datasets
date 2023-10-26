import itertools
import json
import logging
import os
from datetime import datetime

from tqdm import tqdm

from csmed.csmed_cochrane.prepare_dataset import prepare_dataset

if __name__ == "__main__":
    TAR_REPOSITORY_PATH = "../../tar/"

    dataset_names = ["2017-TAR", "2018-TAR", "2019-TAR"]
    review_types = ["Intervention", "Prognosis", "Qualitative", "DTA"]
    dataset_splits = ["Training", "Testing"]

    for dataset_name in dataset_names:
        added_reviews = []
        index_id = 1
        OUTPUT_DATA_PATH = f"../data/external/tar{dataset_name[:4]}/"
        unsuccessful_reviews = []

        if not os.path.exists(OUTPUT_DATA_PATH):
            os.makedirs(OUTPUT_DATA_PATH)

        for review_type, dataset_split in itertools.product(
            review_types, dataset_splits
        ):
            if dataset_name in ["2018-TAR", "2017-TAR"] and review_type != "DTA":
                # 2018-TAR and 2017-TAR only has DTA reviews
                continue

            if dataset_name == "2019-TAR":
                reviews_path = f"{TAR_REPOSITORY_PATH}/{dataset_name}/Task1/{dataset_split}/{review_type}/protocols/"
            elif dataset_name == "2018-TAR":
                reviews_path = f"{TAR_REPOSITORY_PATH}/{dataset_name}/Task1/{dataset_split}/protocols/"
            elif dataset_name == "2017-TAR":
                reviews_path = f"{TAR_REPOSITORY_PATH}/{dataset_name}/{dataset_split.lower()}/extracted_data/"

            if not os.path.exists(reviews_path):
                continue

            print(f"Processing {dataset_name} {review_type} {dataset_split}")
            for review_file in tqdm(os.listdir(reviews_path)):
                if dataset_name == "2017-TAR" and not review_file.endswith(".pids"):
                    continue

                review_id = review_file.split(".")[0]

                tar_dataset = prepare_dataset(
                    review_id=review_id, output_data_path=OUTPUT_DATA_PATH
                )
                if tar_dataset == {}:
                    print(f"No data found for {review_id}")
                    unsuccessful_reviews.append(review_id)
                    continue
                tar_dataset["id"] = index_id
                tar_dataset["used_in"] = [f"{dataset_name} {dataset_split}"]
                tar_dataset["review_type_manual"] = review_type
                if dataset_split == "Training":
                    tar_dataset["dataset_split"] = "train"
                elif dataset_split == "Testing":
                    tar_dataset["dataset_split"] = "test"

                added_reviews.append(tar_dataset)
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
