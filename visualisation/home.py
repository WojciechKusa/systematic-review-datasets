import streamlit as st
import json
import pandas as pd

data_index_file = "data/external/ec2s/data_index.json"

with open(data_index_file) as f:
    data_index = json.load(f)

data_index_df = pd.DataFrame(data_index['data'])

st.dataframe(data_index_df)


