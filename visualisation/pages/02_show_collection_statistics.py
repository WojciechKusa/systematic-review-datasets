import json
import os
import xml.etree.ElementTree as ET

import pandas as pd
import streamlit as st

from visualisation.utils import convert_xml_to_json, get_main_text

DATA_PATH = "data/external/"

folders = [x for x in os.listdir(DATA_PATH) if (os.path.isdir(os.path.join(DATA_PATH, x)) and os.path.isfile(os.path.join(DATA_PATH, x, "data_index.json")))]

st.sidebar.title("Collection selection")
collection = st.sidebar.selectbox("Select a collection", folders)

data_index_file = f"{DATA_PATH}/{collection}/data_index.json"

with open(data_index_file) as f:
    data_index = json.load(f)

data_index_df = pd.DataFrame(data_index["data"])

datasets = data_index_df["review_id"].unique()
st.sidebar.write(f"Collection contains {len(datasets)} datasets")


all_references_df = pd.DataFrame()
for dataset in datasets:
    references_df = pd.read_csv(f"{DATA_PATH}/{collection}/{dataset}/references.csv")
    references_df["dataset"] = dataset
    # df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    # all_references_df = all_references_df.append(references_df)
    #     use concat
    all_references_df = pd.concat([all_references_df, references_df])

st.dataframe(all_references_df)

stats_dict = {}
for reference_type in ["included", "excluded"]:
    stats_dict[reference_type] = {}
    stats_dict[reference_type]["total"] = len(
        all_references_df[all_references_df["reference_type"] == reference_type]
    )
    stats_dict[reference_type]["pubmed_id"] = len(
        all_references_df[
            (all_references_df["reference_type"] == reference_type)
            & (all_references_df["PUBMED_ID"].notnull())
        ]
    )
    stats_dict[reference_type]["doi"] = len(
        all_references_df[
            (all_references_df["reference_type"] == reference_type)
            & (all_references_df["doi"].notnull())
        ]
    )
    stats_dict[reference_type]["pdf"] = len(
        all_references_df[
            (all_references_df["reference_type"] == reference_type)
            & (
                (
                    all_references_df["S2_PDF_LINK"].notnull()
                    | (all_references_df["CORE_PDF_LINK"].notnull())
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
if "exclusion_reason" in all_references_df.columns:
    st.dataframe(
        all_references_df[all_references_df["reference_type"] == "excluded"][
            "exclusion_reason"
        ].value_counts()
    )
else:
    st.write("No exclusion reasons found")

extracted_ft_df = pd.DataFrame()
st.title("Statistics of full-texts")
full_text_files = []
for dataset in datasets:
    full_text_folder = f"{DATA_PATH}/{collection}/{dataset}/pdfs/"
    if os.path.exists(f"{DATA_PATH}/{collection}/{dataset}/extracted_full_texts.csv"):
        extracted_ft_df = pd.concat(
            [
                extracted_ft_df,
                pd.read_csv(f"{DATA_PATH}/{collection}/{dataset}/extracted_full_texts.csv"),
            ],
            ignore_index=True,
        )

    try:
        full_text_files.extend(
            [
                f"{full_text_folder}/{f}"
                for f in os.listdir(full_text_folder)
                if f.endswith(".tei.xml")
            ]
        )
    except FileNotFoundError:
        pass

st.write(f"Found {len(full_text_files)} full-texts")

# only included and excluded references
extracted_ft_df = extracted_ft_df[
    extracted_ft_df["decision"].isin(["included", "excluded"])
]
st.dataframe(extracted_ft_df)
st.write(f"Found {len(extracted_ft_df)} extracted full-texts with decision")
st.write(
    extracted_ft_df["decision"].value_counts() / len(extracted_ft_df) * 100,
    "Percentage of extracted full-texts with decision",
)


avg_full_text_length_characters = 0
avg_full_text_length_words = 0

for file_path in full_text_files:
    title, abstract, main_text = get_main_text(file_path)

    full_text = title + abstract + main_text

    avg_full_text_length_characters += len(full_text)
    avg_full_text_length_words += len(full_text.split())
    tree = ET.parse(file_path)
    root = tree.getroot()
    json_data = convert_xml_to_json(file_path)

avg_full_text_length_characters = int(
    avg_full_text_length_characters / len(full_text_files)
)
avg_full_text_length_words = int(avg_full_text_length_words / len(full_text_files))
st.write(
    f"Average full-text length: {avg_full_text_length_characters} characters,"
    f" {avg_full_text_length_words} words, and {avg_full_text_length_words * 1.5} tokens"
)
