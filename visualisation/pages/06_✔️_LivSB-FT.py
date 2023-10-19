import os
import zipfile

import pandas as pd
import streamlit as st

data_path = os.path.join(os.path.dirname(__file__), "../../data/CSMeD/CSMeD-FT/")

if not os.path.exists(data_path):
    st.write("Data not available. Unzipping data...")

    with zipfile.ZipFile(
        os.path.join(os.path.dirname(__file__), "../../data/CSMeD/CSMeD-FT.zip"), "r"
    ) as zip_ref:
        zip_ref.extractall(os.path.join(os.path.dirname(__file__), "../../data/CSMeD/"))

    st.write("Done unzipping data.")

datasets = [x for x in os.listdir(data_path) if x.endswith(".csv")]

dfs = {}
for dataset in datasets:
    dfs[dataset] = pd.read_csv(f"{data_path}/{dataset}")

dataset = st.sidebar.selectbox(
    "Select a dataset",
    datasets,
    index=0,
)

data2 = st.sidebar.selectbox(
    "Select another dataset",
    datasets,
    index=1,
)

# measure overlap between two datasets in review_ids
review_ids1 = set(dfs[dataset]["review_id"].unique())
review_ids2 = set(dfs[data2]["review_id"].unique())

st.write(f"Number of reviews in {dataset}: ", len(review_ids1))
st.write(f"Number of reviews in {data2}: ", len(review_ids2))

st.write(
    f"Number of reviews in both datasets: ", len(review_ids1.intersection(review_ids2))
)
st.write(
    f"Number of reviews in {dataset} but not in {data2}: ",
    len(review_ids1.difference(review_ids2)),
)
st.write(
    f"Number of reviews in {data2} but not in {dataset}: ",
    len(review_ids2.difference(review_ids1)),
)
