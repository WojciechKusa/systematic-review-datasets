from collections import Counter

import numpy as np

from csmed.datasets.datasets.csmed_basic.csmed_basic import CSMeDBasic
from csmed.datasets.datasets.csmed_cochrane.csmed_cochrane import CSMeDCochrane


def print_statistics_cochrane(dataset) -> None:
    systematic_review_names = dataset.keys()

    # average length of abstracts
    abstract_lengths = [
        len(dataset[key]["dataset_details"]["abstract"].split())
        for key in systematic_review_names
    ]
    avg_abstract_length = sum(abstract_lengths) / len(abstract_lengths)
    print(f"average length of abstracts: {avg_abstract_length}")

    # average length of titles
    title_lengths = [
        len(dataset[key]["dataset_details"]["title"].split())
        for key in systematic_review_names
    ]
    avg_title_length = sum(title_lengths) / len(title_lengths)
    print(f"average length of titles: {avg_title_length}")

    # average length of search queries
    query_lengths = []
    for key in systematic_review_names:
        try:
            review_query = list(
                dataset[key]["dataset_details"]["search_strategy"].values()
            )[0]
        except IndexError:
            review_query = "no search query"
        query_lengths.append(len(review_query.split()))
    avg_query_length = sum(query_lengths) / len(query_lengths)
    print(f"average length of queries: {avg_query_length}")

    # average length of documents criteria
    criteria_lengths = [
        " ".join(
            [
                f"{k}: {v}"
                for k, v in dataset[key]["dataset_details"]["criteria"].items()
            ]
        )
        for key in systematic_review_names
    ]
    criteria_lengths = [len(x.split()) for x in criteria_lengths]
    avg_criteria_length = sum(criteria_lengths) / len(criteria_lengths)
    print(f"average length of criteria: {avg_criteria_length}")

    # calculate medians and percentiles
    abstract_lengths = np.array(abstract_lengths)
    title_lengths = np.array(title_lengths)
    query_lengths = np.array(query_lengths)
    criteria_lengths = np.array(criteria_lengths)

    print(f"median length of abstracts: {np.median(abstract_lengths)}")
    print(f"median length of titles: {np.median(title_lengths)}")
    print(f"median length of queries: {np.median(query_lengths)}")
    print(f"median length of criteria: {np.median(criteria_lengths)}")

    print(f"25th percentile length of abstracts: {np.percentile(abstract_lengths, 25)}")
    print(f"25th percentile length of titles: {np.percentile(title_lengths, 25)}")
    print(f"25th percentile length of queries: {np.percentile(query_lengths, 25)}")
    print(f"25th percentile length of criteria: {np.percentile(criteria_lengths, 25)}")

    print(f"75th percentile length of abstracts: {np.percentile(abstract_lengths, 75)}")
    print(f"75th percentile length of titles: {np.percentile(title_lengths, 75)}")
    print(f"75th percentile length of queries: {np.percentile(query_lengths, 75)}")
    print(f"75th percentile length of criteria: {np.percentile(criteria_lengths, 75)}")

    print(f"90th percentile length of abstracts: {np.percentile(abstract_lengths, 90)}")
    print(f"90th percentile length of titles: {np.percentile(title_lengths, 90)}")
    print(f"90th percentile length of queries: {np.percentile(query_lengths, 90)}")
    print(f"90th percentile length of criteria: {np.percentile(criteria_lengths, 90)}")

    print(f"99th percentile length of abstracts: {np.percentile(abstract_lengths, 99)}")
    print(f"99th percentile length of titles: {np.percentile(title_lengths, 99)}")
    print(f"99th percentile length of queries: {np.percentile(query_lengths, 99)}")
    print(f"99th percentile length of criteria: {np.percentile(criteria_lengths, 99)}")


def print_statistics(dataset) -> None:
    systematic_review_names = dataset.keys()

    sizes = [dataset[key]["data"].num_rows["train"] for key in systematic_review_names]
    n_relevant = [
        Counter(x["labels"][0] for x in dataset[key]["data"]["train"])[RELEVANT_LABEL]
        for key in systematic_review_names
    ]

    print(sizes)
    print(f"total: {sum(sizes)}")
    print(n_relevant)
    print(f"total relevant: {sum(n_relevant)}")

    avg_relevant_per_dataset = sum([x / y for x, y in zip(n_relevant, sizes)]) / len(
        n_relevant
    )
    print(f"average relevant documents per dataset: {avg_relevant_per_dataset}")

    avg_docs_per_dataset = sum(sizes) / len(sizes)
    print(f"average number of documents per dataset: {avg_docs_per_dataset}")

    total_avg_relevant_docs = sum(n_relevant) / sum(sizes)
    print(f"total average number of relevant documents: {total_avg_relevant_docs}")

    lengths = [
        len(doc["text"].split())
        for key in systematic_review_names
        for doc in dataset[key]["data"]["train"]
    ]
    avg_doc_length = sum(lengths) / len(lengths)
    print(f"average length of documents: {avg_doc_length}")


if __name__ == "__main__":
    RELEVANT_LABEL = "1"

    csm = CSMeDBasic()
    csmed_dataset = csm.load_dataset(
        base_path="../../csmed/datasets/datasets/csmed_basic/"
    )
    print_statistics(csmed_dataset["TRAIN"])

    csmed_cochrane = CSMeDCochrane().load_dataset(
        base_path="../../csmed/datasets/datasets/csmed_cochrane/"
    )

    print_statistics_cochrane(csmed_cochrane["TRAIN"])
    print_statistics_cochrane(csmed_cochrane["EVAL"])

    all_datasets = csmed_cochrane["TRAIN"]
    all_datasets.update(csmed_cochrane["EVAL"])

    print_statistics_cochrane(all_datasets)
