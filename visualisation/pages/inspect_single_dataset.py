import streamlit as st
import json
import pandas as pd

DATA_PATH = "data/external/ec2s"

data_index_file = f"{DATA_PATH}/data_index.json"

with open(data_index_file) as f:
    data_index = json.load(f)

data_index_df = pd.DataFrame(data_index['data'])

st.sidebar.title("Select a dataset")
dataset = st.sidebar.selectbox(
    "Select a dataset",
    data_index_df['review_id'].unique()
)

references_df = pd.read_csv(f"{DATA_PATH}/{dataset}/references.csv")

st.dataframe(references_df)


stats_dict = {}
for reference_type in ['included', 'excluded']:
    stats_dict[reference_type] = {}
    stats_dict[reference_type]['total'] = len(references_df[references_df['reference_type'] == reference_type])
    stats_dict[reference_type]['pubmed_id'] = len(references_df[(references_df['reference_type'] == reference_type) & (references_df['PUBMED_ID'].notnull())])
    stats_dict[reference_type]['doi'] = len(references_df[(references_df['reference_type'] == reference_type) & (references_df['doi'].notnull())])
    stats_dict[reference_type]['pdf'] = len(references_df[(references_df['reference_type'] == reference_type) & (references_df['PDF_LINK'].notnull())])

stats_df = pd.DataFrame(stats_dict).T
stats_df['pubmed_id'] = stats_df['pubmed_id'] / stats_df['total']
stats_df['doi'] = stats_df['doi'] / stats_df['total']
stats_df['pdf'] = stats_df['pdf'] / stats_df['total']

stats_df.columns = ['Total', '% with PubMed ID', '% with DOI', '% with PDF']
st.dataframe(stats_df)

st.title("Reasons for exclusion")
st.dataframe(references_df[references_df['reference_type'] == 'excluded']['exclusion_reason'].value_counts())