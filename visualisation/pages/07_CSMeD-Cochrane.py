from collections import Counter

import streamlit as st

from csmed.datasets.datasets.csmed_basic.csmed_basic import CSMeDBasic
from csmed.datasets.datasets.csmed_cochrane.csmed_cochrane import CSMeDCochrane


def count_statistics(dataset):
    systematic_review_names = dataset.keys()

    sizes = [dataset[key]["data"].num_rows["train"] for key in systematic_review_names]
    n_relevant = [
        Counter(x["labels"][0] for x in dataset[key]["data"]["train"])[RELEVANT_LABEL]
        for key in systematic_review_names
    ]

    avg_relevant_per_dataset = sum([x / y for x, y in zip(n_relevant, sizes)]) / len(
        n_relevant
    )

    avg_docs_per_dataset = sum(sizes) / len(sizes)

    total_avg_relevant_docs = sum(n_relevant) / sum(sizes)

    lengths = [
        len(doc["text"].split())
        for key in systematic_review_names
        for doc in dataset[key]["data"]["train"]
    ]

    stats = {
        "avg_relevant_per_dataset": avg_relevant_per_dataset,
        "avg_docs_per_dataset": avg_docs_per_dataset,
        "total_avg_relevant_docs": total_avg_relevant_docs,
        "avg_doc_length": sum(lengths) / len(lengths),
        "total_size": sum(sizes),
        "total_relevant": sum(n_relevant),
    }
    return stats


RELEVANT_LABEL = "1"

csmed_dataset = CSMeDBasic().load_dataset(
    base_path="csmed/datasets/datasets/csmed_basic/"
)

stats = count_statistics(csmed_dataset["TRAIN"])
st.write(stats)

csmed_cochrane = CSMeDCochrane().load_dataset(
    base_path="csmed/datasets/datasets/csmed_cochrane/"
)

stats = count_statistics(csmed_cochrane["TRAIN"])
st.write(stats)

st.write("EVAL")
stats = count_statistics(csmed_cochrane["EVAL"])
st.write(stats)
