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
            "collection": collection,
            "#datasets": len(collections[collection]),
            "#included FT": sum(
                [
                    dataset["n_references_included"]
                    for dataset in collections[collection]
                ]
            ),
            "avg #included FT": sum(
                [
                    dataset["n_references_included"]
                    for dataset in collections[collection]
                ]
            )
            / len(collections[collection]),
            # "#excluded": sum([len(dataset["excluded"]) for dataset in collections[collection]]),
        }
    )
stats_df = pd.DataFrame(stats)
st.dataframe(stats_df)

with st.expander("Review overlap"):
    review_counts = {}
    for collection in collections:
        for review in collections[collection]:
            if review["review_id"] not in review_counts:
                review_counts[review["review_id"]] = []
            review_counts[review["review_id"]].append(collection)

    for review in review_counts:
        if len(review_counts[review]) > 1:
            st.write(review, review_counts[review])
