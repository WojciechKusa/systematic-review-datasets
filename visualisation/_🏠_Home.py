import json
import os

import pandas as pd
import streamlit as st

st.title(
    "This is the visualisation for the LivSB meta-dataset. Use the sidebar for navigation."
)

st.write(
    "Some of the data is not available in the repository, but can be downloaded from the original source."
    "These pages will not work on the deployed version of the app, but can be run locally."
)

st.write(
    "Pages available in the deployed version of the app are marked with a ✔️ emoji."
)

data_index_file = os.path.join(
    os.path.dirname(__file__), "../data/external/ec2s/data_index.json"
)

with open(data_index_file) as f:
    data_index = json.load(f)

data_index_df = pd.DataFrame(data_index["data"])

st.dataframe(data_index_df)
