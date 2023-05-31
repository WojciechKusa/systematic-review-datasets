import json
import os
import xml.etree.ElementTree as ET

import pandas as pd
import streamlit as st

from visualisation.utils import convert_xml_to_json, get_main_text


DATA_PATH = "data/external/"

folders = [x for x in os.listdir(DATA_PATH) if (os.path.isdir(os.path.join(DATA_PATH, x)) and os.path.isfile(os.path.join(DATA_PATH, x, "data_index.json")))]

st.sidebar.title("Review selection")
collection = st.sidebar.selectbox("Select a collection", folders)

DATA_PATH = f"{DATA_PATH}/{collection}"

data_index_file = f"{DATA_PATH}/data_index.json"

with open(data_index_file) as f:
    data_index = json.load(f)

data_index_df = pd.DataFrame(data_index["data"])

dataset = st.sidebar.selectbox("Select a review", data_index_df["review_id"].unique())

references_df = pd.read_csv(f"{DATA_PATH}/{dataset}/references.csv")

st.dataframe(references_df)

stats_dict = {}
for reference_type in ["included", "excluded"]:
    stats_dict[reference_type] = {}
    stats_dict[reference_type]["total"] = len(
        references_df[references_df["reference_type"] == reference_type]
    )
    stats_dict[reference_type]["pubmed_id"] = len(
        references_df[
            (references_df["reference_type"] == reference_type)
            & (references_df["PUBMED_ID"].notnull())
        ]
    )
    stats_dict[reference_type]["doi"] = len(
        references_df[
            (references_df["reference_type"] == reference_type)
            & (references_df["doi"].notnull())
        ]
    )
    stats_dict[reference_type]["pdf"] = len(
        references_df[
            (references_df["reference_type"] == reference_type)
            & (
                (
                    references_df["S2_PDF_LINK"].notnull()
                    | (references_df["CORE_PDF_LINK"].notnull())
                )
            )
        ]
    )

stats_df = pd.DataFrame(stats_dict).T
stats_df["pubmed_id"] = stats_df["pubmed_id"] / stats_df["total"]
stats_df["doi"] = stats_df["doi"] / stats_df["total"]
stats_df["pdf"] = stats_df["pdf"] / stats_df["total"]

stats_df.columns = ["Total", "% with PubMed ID", "% with DOI", "% with PDF"]
st.dataframe(stats_df)

st.title("Reasons for exclusion")
if "exclusion_reason" in references_df.columns:
    st.dataframe(
        references_df[references_df["reference_type"] == "excluded"][
            "exclusion_reason"
        ].value_counts()
    )
else:
    st.write("No exclusion reasons found")

st.title("Statistics of full-texts")
full_text_folder = f"{DATA_PATH}/{dataset}/pdfs/"
full_text_files = [f for f in os.listdir(full_text_folder) if f.endswith(".tei.xml")]
st.write(f"Found {len(full_text_files)} full-texts")

avg_full_text_length_characters = 0
avg_full_text_length_words = 0

for file in full_text_files:
    title, abstract, main_text = get_main_text(f"{full_text_folder}/{file}")

    full_text = title + abstract + main_text

    avg_full_text_length_characters += len(full_text)
    avg_full_text_length_words += len(full_text.split())
    tree = ET.parse(f"{full_text_folder}/{file}")
    root = tree.getroot()
    json_data = convert_xml_to_json(f"{full_text_folder}/{file}")

avg_full_text_length_characters = int(
    avg_full_text_length_characters / len(full_text_files)
)
avg_full_text_length_words = int(avg_full_text_length_words / len(full_text_files))
st.write(
    f"Average full-text length: {avg_full_text_length_characters} characters,"
    f" {avg_full_text_length_words} words, and {avg_full_text_length_words*1.4} tokens"
)

file = st.selectbox("Select a full-text", full_text_files)
title, abstract, main_text = get_main_text(f"{full_text_folder}/{file}")

st.write(f"**Title**: {title}")
st.write(f"**Abstract**: {abstract}")
st.write(f"**Main text**: {main_text}")
