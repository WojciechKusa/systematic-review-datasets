from ranx import Qrels, Run
from ranx import compare
from retriv import SparseRetriever, DenseRetriever
import time

from csmed.datasets.datasets.csmed_cochrane.csmed_cochrane import CSMeDCochrane


def evaluate_runs(
    runs: dict[str, dict[str, dict[str, float]]],
    qrels_dict: dict[str, dict[str, int]],
    outfile: str = "report.tex",
) -> None:
    qrels = Qrels(qrels_dict)
    runs_list = [Run(_run, name=_run_name) for _run_name, _run in runs.items()]

    report = compare(
        qrels=qrels,
        runs=runs_list,
        metrics=[
            "map@100",
            "mrr@100",
            "ndcg@10",
            "map",
            "recall@10",
            "recall@50",
            "recall@100",
            "recall",
            "precision",
            "precision@10",
            "precision@50",
            "precision@100",
            "precision",
            "f1",
            "f1@10",
            "f1@50",
            "f1@100",
        ],
        max_p=0.05,
    )

    print(report)
    with open(outfile, "w", encoding="utf_8") as f:
        f.write(report.to_latex())


if __name__ == "__main__":
    csmed_dataset = CSMeDCochrane().load_dataset(
        base_path="../../csmed/datasets/datasets/csmed_cochrane/"
    )

    eval_reviews = csmed_dataset["EVAL"]

    qrels_dict = {}

    runs = {
        "bm25_title": {},
        "bm25_abstract": {},
        "bm25_criteria": {},
        "bm25_query": {},
    }

    for index_i, (review_name, review_data) in enumerate(eval_reviews.items(), start=1):
        collection = [
            {"id": doc["document_id"], "text": doc["text"]}
            for doc in review_data["data"]["train"]
        ]
        qrels = {
            doc["document_id"]: int(doc["labels"][0])
            for doc in review_data["data"]["train"]
        }

        qrels_dict[review_name] = qrels

        sr = SparseRetriever(
            index_name="new-index",
            model="bm25",
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
        sr.index(collection)


        review_title = review_data["dataset_details"]["title"]
        review_abstract = review_data["dataset_details"]["abstract"]
        review_criteria = " ".join(review_data["dataset_details"]["criteria"])
        try:
            review_query = list(
                review_data["dataset_details"]["search_strategy"].values()
            )[0]
        except IndexError:
            review_query = "no search query"
        collection_size = len(review_data["data"]["train"])

        runs["bm25_criteria"][review_name] = sr.search(
            query=review_criteria, cutoff=collection_size, return_docs=False
        )
        runs["bm25_query"][review_name] = sr.search(
            query=review_query, cutoff=collection_size, return_docs=False
        )
        runs["bm25_title"][review_name] = sr.search(
            query=review_title, cutoff=collection_size, return_docs=False
        )
        runs["bm25_abstract"][review_name] = sr.search(
            query=review_abstract, cutoff=collection_size, return_docs=False
        )
        if index_i % 10 == 0:
            evaluate_runs(runs=runs, qrels_dict=qrels_dict)

    evaluate_runs(runs=runs, qrels_dict=qrels_dict)
