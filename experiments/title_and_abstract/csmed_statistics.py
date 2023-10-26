from collections import Counter

from csmed.datasets.datasets.csmed_basic.csmed_basic import CSMeDBasic

if __name__ == "__main__":
    RELEVANT_LABEL = "1"

    csm = CSMeDBasic()
    csmed_dataset = csm.load_dataset(
        base_path="../../csmed/datasets/datasets/csmed_basic/"
    )

    systematic_review_names = csmed_dataset.keys()

    sizes = [csmed_dataset[key].num_rows["train"] for key in systematic_review_names]
    n_relevant = [
        Counter(x["labels"][0] for x in csmed_dataset[key]["train"])[RELEVANT_LABEL]
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
        for doc in csmed_dataset[key]["train"]
    ]
    avg_doc_length = sum(lengths) / len(lengths)
    print(f"average length of documents: {avg_doc_length}")
