import json

from ranx import Qrels, Run, compare


def n_precision_at_recall(
    run: dict[str, dict[str, float]],
    qrels_dict: dict[str, dict[str, int]],
    recall_level: float = 0.95,
) -> float:
    """calculates normalised Precision achieved at a given Recall level.
    Recall level is achieved when a greater or equal number of docs obtains that recall.
        "normalisedPrecision": r"nPrecision@r\% &= \frac{TP \cdot TN}{\mathcal{E} \cdot (TP + FP)} ",

    """
    n_precision_values = []

    for query_id, rankings in run.items():
        tp, fp, tn, fn = 0, 0, 0, 0

        relevant_docs = sum(qrels_dict[query_id].values())
        total_docs = len(qrels_dict[query_id])
        non_relevant_docs = total_docs - relevant_docs

        for rank, (doc_id, _) in enumerate(rankings.items()):
            if qrels_dict[query_id].get(doc_id, 0) == 1:  # Document is relevant
                tp += 1
            else:
                fp += 1

            recall = tp / (tp + fn + relevant_docs - tp)

            if recall >= recall_level:
                tn = non_relevant_docs - fp
                fn = relevant_docs - tp
                # E = FP + TN
                # n_precision = (TP * TN) / (E * (TP + FP)) if (E * (TP + FP)) != 0 else 0
                n_precision = (tp * tn) / ((fp + tn) * (tp + fp))
                n_precision_values.append(n_precision)
                break

    return (
        sum(n_precision_values) / len(n_precision_values) if n_precision_values else 0
    )


def tnr_at_recall(
    run: dict[str, dict[str, float]],
    qrels_dict: dict[str, dict[str, int]],
    recall_level: float = 0.95,
) -> float:
    """calculates the True Negative Rate achieved at a given Recall level.
    Recall level is achieved when a greater or equal number of docs obtains that recall.
    """
    tnr_values = []

    for query_id, rankings in run.items():
        tp, fp, tn, fn = 0, 0, 0, 0

        relevant_docs = sum(qrels_dict[query_id].values())
        total_docs = len(qrels_dict[query_id])
        non_relevant_docs = total_docs - relevant_docs

        for rank, (doc_id, _) in enumerate(rankings.items()):
            if qrels_dict[query_id].get(doc_id, 0) == 1:  # Document is relevant
                tp += 1
            else:
                fp += 1

            recall = tp / (tp + fn + relevant_docs - tp)

            if recall >= recall_level:
                tn = non_relevant_docs - fp
                fn = relevant_docs - tp
                tnr = tn / (tn + fp)
                tnr_values.append(tnr)
                break

    return sum(tnr_values) / len(tnr_values) if tnr_values else 0


def find_last_relevant(
    run: dict[str, dict[str, float]], qrels_dict: dict[str, dict[str, int]]
) -> float:
    """Find percentage of the run where the last relevant item was found."""
    percentages = []

    for query, docs in run.items():
        # Check if the query exists in the qrels_dict
        if query in qrels_dict:
            # List of relevant document ids for the given query
            relevant_docs = [
                doc_id for doc_id, rel in qrels_dict[query].items() if rel > 0
            ]

            # Find the position of the last relevant document
            last_relevant_position = None
            for doc_id, _ in reversed(docs.items()):
                if doc_id in relevant_docs:
                    last_relevant_position = (
                        list(docs.keys()).index(doc_id) + 1
                    )  # Adding 1 as indexing starts from 0
                    break

            # If a relevant document is found in the run
            if last_relevant_position is not None:
                percentages.append(last_relevant_position / len(docs) * 100)

    # Return the average of the percentages
    return sum(percentages) / len(percentages) if percentages else 0.0


def evaluate_runs(
    runs: dict[str, dict[str, dict[str, float]]],
    qrels_dict: dict[str, dict[str, int]],
    outfile: str = "report.tex",
    other_measures_file: str = "other_measures.json",
    outfile_path: str = "../../reports/title_and_abstract/",
) -> None:
    bonferroni_correction = 0.05 / len(runs)

    qrels = Qrels(qrels_dict)
    runs_list = [Run(_run, name=_run_name) for _run_name, _run in runs.items()]

    report = compare(
        qrels=qrels,
        runs=runs_list,
        metrics=[
            "ndcg",
            "ndcg@5",
            "ndcg@10",
            "ndcg@100",
            "map",
            "map@10",
            "map@100",
            "recall",
            "recall@10",
            "recall@50",
            "recall@100",
            "precision",
            "precision@10",
            "precision@50",
            "precision@100",
            "f1",
            "f1@10",
            "f1@50",
            "f1@100",
            "r-precision",
            "mrr@100",
        ],
        max_p=bonferroni_correction,
    )

    print(report)
    with open(f"{outfile_path}/{outfile}", "w") as f:
        f.write(report.to_latex())

    other_measures = {
        "last_rel": {},
        "tnr@95": {},
        "tnr@90": {},
        "n_precision@95": {},
        "n_precision@80": {},
        "dataset_size": {},
    }
    for run_name, run in runs.items():
        other_measures["last_rel"][run_name] = find_last_relevant(
            run=run, qrels_dict=qrels_dict
        )
        other_measures["tnr@95"][run_name] = tnr_at_recall(
            run=run, qrels_dict=qrels_dict, recall_level=0.95
        )
        other_measures["tnr@90"][run_name] = tnr_at_recall(
            run=run, qrels_dict=qrels_dict, recall_level=0.90
        )
        other_measures["n_precision@95"][run_name] = n_precision_at_recall(
            run=run, qrels_dict=qrels_dict, recall_level=0.95
        )
        other_measures["n_precision@80"][run_name] = n_precision_at_recall(
            run=run, qrels_dict=qrels_dict, recall_level=0.80
        )
        # other_measures["dataset_size"][run_name] = []
    print(other_measures)
    with open(f"{outfile_path}/{other_measures_file}", "w") as f:
        json.dump(other_measures, f, indent=4)
