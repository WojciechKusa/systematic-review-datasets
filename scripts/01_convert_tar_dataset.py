import itertools
import json
import os
from datetime import datetime

from tqdm import tqdm

from ec2s.ec2s.prepare_dataset import prepare_dataset

if __name__ == "__main__":
    TAR_REPOSITORY_PATH = "../../tar/"
    OUTPUT_DATA_PATH = "../data/external/tar-systematic_reviews/"

    if not os.path.exists(OUTPUT_DATA_PATH):
        os.makedirs(OUTPUT_DATA_PATH)

    dataset_names = ["2017-TAR", "2018-TAR", "2019-TAR"]
    review_types = ["Intervention", "Prognosis", "Qualitative", "DTA"]
    dataset_splits = ["Training", "Testing"]

    added_reviews = []

    index_id = 1

    for dataset_name in dataset_names:
        for review_type, dataset_split in itertools.product(review_types, dataset_splits):
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

            for review_file in tqdm(os.listdir(reviews_path)):
                if dataset_name == "2017-TAR" and not review_file.endswith(".pids"):
                    continue

                review_id = review_file.split(".")[0]

                tar_dataset = prepare_dataset(
                    review_id=review_id, output_data_path=OUTPUT_DATA_PATH
                )
                tar_dataset["id"] = index_id
                tar_dataset["used_in"] = [f"{dataset_name} {dataset_split}"]
                tar_dataset["review_type"] = review_type
                tar_dataset["dataset_split"] = dataset_split

                added_reviews.append(tar_dataset)
                index_id += 1

                final_data = {
                    "data_format_version": "1.0.0",
                    "updated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "data": added_reviews,
                }
                with open(f"{OUTPUT_DATA_PATH}/data_index.json", "w") as f:
                    json.dump(final_data, f, indent=2)
