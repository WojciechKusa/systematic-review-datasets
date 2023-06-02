import datasets
import streamlit as st


available_datasets = [
    "cohen",
    "swift"
]

dataset_name = st.sidebar.selectbox(
    "Select a dataset",
    available_datasets
)

x = datasets.load_dataset(dataset_name)

# show all versions of the dataset
versions = datasets.list_versions(dataset_name)


st.title("cohen")
st.write(x)