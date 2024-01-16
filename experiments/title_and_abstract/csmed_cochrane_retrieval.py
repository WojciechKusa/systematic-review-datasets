import os.path
import pickle
import random
from typing import Any

import numpy as np
import torch
from numba import njit
from retriv import SparseRetriever, DenseRetriever

from csmed.datasets.datasets.csmed_cochrane.csmed_cochrane import CSMeDCochrane
from experiments.title_and_abstract.measures import evaluate_runs

# Constants
SEED = 42
USE_GPU = True

# Initialize seed for reproducibility
torch.manual_seed(SEED)
random.seed(SEED)
np.random.seed(SEED)


@njit
def set_seed(value):
    np.random.seed(value)


def load_dataset():
    return CSMeDCochrane().load_dataset(
        base_path="../../csmed/datasets/datasets/csmed_cochrane/"
    )


def create_retrievers(configs, collection):
    retrievers = {}
    for name, conf in configs.items():
        if conf["type"] == "sparse":
            retrievers[name] = SparseRetriever(
                index_name=f"new-index{conf['model']}",
                model=conf["model"],
                min_df=1,
                tokenizer="whitespace",
                stemmer="english",
                stopwords="english",
                do_lowercasing=True,
                do_ampersand_normalization=True,
                do_special_chars_normalization=True,
                do_acronyms_normalization=True,
                do_punctuation_removal=True,
            )
            retrievers[name].index(collection)
        elif conf["type"] == "dense":
            retrievers[name] = DenseRetriever(
                index_name=f"dense_{conf['model']}",
                model=conf["model"],
                normalize=True,
                max_length=conf["max_length"],
                use_ann=True,
            )
            retrievers[name].index(collection, use_gpu=USE_GPU, batch_size=128)
    return retrievers


def extract_review_details(review_data):
    review_details = {
        "title": review_data["dataset_details"]["title"],
        "abstract": review_data["dataset_details"]["abstract"],
        "criteria": " ".join(
            [f"{k}: {v}" for k, v in review_data["dataset_details"]["criteria"].items()]
        ),
    }

    # Handling the case where the search strategy might not be present
    try:
        review_details["query"] = list(
            review_data["dataset_details"]["search_strategy"].values()
        )[0]
    except IndexError:
        review_details["query"] = "no search query"

    return review_details


def prepare_data(review_data):
    # Extracting documents for indexing
    collection = [
        {"id": doc["document_id"], "text": doc["text"]}
        for doc in review_data["data"]["train"]
    ]

    # Extracting relevance judgments (qrels)
    qrels = {
        doc["document_id"]: int(doc["labels"][0])
        for doc in review_data["data"]["train"]
    }

    return collection, qrels


def initialise_runs(retriever_configs: dict[str, dict[str, str]]) -> dict[str, Any]:
    runs = {}
    for retriever_name in retriever_configs.keys():
        for query_type in ["abstract", "title", "query", "criteria"]:
            run_key = f"{retriever_name}_{query_type}"
            runs[run_key] = {}
    return runs


def process_review(
    retrievers,
    review_data,
    runs: dict[str, dict[str, dict[str, float]]],
    collection_size: int,
    review_name: str,
):
    review_details = extract_review_details(review_data)

    for model_name, retriever in retrievers.items():
        for slr_protocol_key, slr_protocol_value in review_details.items():
            run_key = f"{model_name}_{slr_protocol_key}"
            # print(slr_protocol_value)
            # print(run_key)
            runs[run_key][review_name] = retriever.search(
                query=slr_protocol_value, cutoff=collection_size, return_docs=False
            )


if __name__ == "__main__":
    set_seed(SEED)

    dataset = load_dataset()
    eval_reviews = dataset["EVAL"]

    outfile_path = "../../reports/title_and_abstract/"
    if not os.path.exists(outfile_path):
        os.makedirs(outfile_path)

    retriever_configs = {
        "bm25": {
            "type": "sparse",
            "model": "bm25",
        },
        "tf-idf": {
            "type": "sparse",
            "model": "tf-idf",
        },
        # "MiniLM-128": {
        #     "type": "dense",
        #     "model": "sentence-transformers/all-MiniLM-L6-v2",
        #     "max_length": 128,
        # },
        "MiniLM-256": {
            "type": "dense",
            "model": "sentence-transformers/all-MiniLM-L6-v2",
            "max_length": 256,
        },
        # "qa-MiniLM-512": {
        #     "type": "dense",
        #     "model": "sentence-transformers/multi-qa-MiniLM-L6-cos-v1",
        #     "max_length": 512,
        # },
        "mpnet": {
            "type": "dense",
            "model": "sentence-transformers/all-mpnet-base-v2",
            "max_length": 512,
        },
        # "nli-mpnet": {
        #     "type": "dense",
        #     "model": "sentence-transformers/nli-mpnet-base-v2",
        #     "max_length": 512,
        # },
        # "biobert-nli": {
        #     "type": "dense",
        #     "model": "pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb",
        #     "max_length": 512,
        # },
        "S-BioBert": {
            "type": "dense",
            "model": "pritamdeka/S-BioBert-snli-multinli-stsb",
            "max_length": 512,
        },
        "pubmedbert": {
            "type": "dense",
            "model": "pritamdeka/S-PubMedBert-MS-MARCO",
            "max_length": 512,
        },
        # "roberta": {
        #     "type": "dense",
        #     "model": "sentence-transformers/stsb-roberta-base-v2",
        #     "max_length": 512,
        # },
    }

    dataset = load_dataset()
    eval_reviews = dataset["EVAL"]
    qrels_dict = {}

    runs = initialise_runs(retriever_configs)

    print(f"Processing {len(eval_reviews)} reviews with {len(runs)} runs")
    for index, (review_name, review_data) in enumerate(eval_reviews.items(), start=1):
        collection, qrels = prepare_data(review_data)
        qrels_dict[review_name] = qrels

        retrievers = create_retrievers(retriever_configs, collection=collection)
        process_review(
            retrievers,
            review_data,
            runs,
            collection_size=len(review_data["data"]["train"]),
            review_name=review_name,
        )
        if index % 10 == 0:
            evaluate_runs(runs=runs, qrels_dict=qrels_dict)

            with open(f"{outfile_path}/runs.pkl", "wb") as f:
                pickle.dump(runs, f)

    evaluate_runs(runs=runs, qrels_dict=qrels_dict)
    with open(f"{outfile_path}/runs.pkl", "wb") as f:
        pickle.dump(runs, f)

