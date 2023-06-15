import json
import os
import xml.etree.ElementTree as ET

import pandas as pd
import streamlit as st

from visualisation.utils import convert_xml_to_json, get_main_text

DATA_PATH = "data/external/"

folders = [
    x
    for x in os.listdir(DATA_PATH)
    if (
        os.path.isdir(os.path.join(DATA_PATH, x))
        and os.path.isfile(os.path.join(DATA_PATH, x, "data_index.json"))
    )
]

collections = {}
for folder in folders:
    data_index_file = f"{DATA_PATH}/{folder}/data_index.json"
    with open(data_index_file) as f:
        data_index = json.load(f)
    collections[folder] = data_index["data"]

# write a table with all the dataset statistics: #reviews, average #included publications
st.title("Dataset statistics")
st.write("This page shows the statistics for all the datasets in the corpus.")
stats = []
for collection in collections:
    stats.append(
        {
            "Collection name": collection,
            "# reviews": len(collections[collection]),
            "# unique reviews": len(
                list(set([x["review_id"] for x in collections[collection]]))
            ),
            "# reviews in train": sum(
                [1 for x in collections[collection] if x["dataset_split"] == "train"]
            ),
            "# reviews in test": sum(
                [1 for x in collections[collection] if x["dataset_split"] == "test"]
            ),
            "# included FT": sum(
                [
                    dataset["unique_references_included"]
                    for dataset in collections[collection]
                ]
            ),
            "Avg. #included FT": int(
                sum(
                    [
                        dataset["unique_references_included"]
                        for dataset in collections[collection]
                    ]
                )
                / len(collections[collection])
            ),
            # "# excluded FT": sum(
            #     [
            #         dataset["unique_references_excluded"]
            #         for dataset in collections[collection]
            #     ]
            # ),
            # "avg #excluded FT": int(sum(
            #     [
            #         dataset["unique_references_excluded"]
            #         for dataset in collections[collection]
            #     ]
            # )
            # / len(collections[collection])),
        }
    )
stats_df = pd.DataFrame(stats)
st.dataframe(stats_df)

with st.expander("Review overlap"):
    review_counts = {}
    for collection in sorted(collections):
        for review in collections[collection]:
            if review["review_id"] not in review_counts:
                review_counts[review["review_id"]] = []
            review_counts[review["review_id"]].append(
                f"{collection} ({review['dataset_split']})"
            )

    overlapping_reviews = []
    for review in review_counts:
        if len(review_counts[review]) > 1:
            st.write(review, review_counts[review])
            overlapping_reviews.append(review)

overlapping_reviews_df = []
for review in overlapping_reviews:

    overlapping_reviews_df.append(
        {
            "review_id": review,
            "first_collection": review_counts[review][0],
            "other_collections": ", ".join(review_counts[review][1:]),
        }
    )
overlapping_reviews_df = pd.DataFrame(overlapping_reviews_df)
st.dataframe(overlapping_reviews_df)

st.write(f"Total number of overlapping reviews: {len(overlapping_reviews)}")

st.write(f"Total number of unique reviews: {len(review_counts)}")
