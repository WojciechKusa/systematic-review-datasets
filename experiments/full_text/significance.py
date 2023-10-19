import json
import os

import numpy as np
import pandas as pd
from datasets import load_dataset
from sklearn.dummy import DummyClassifier
from statsmodels.stats.contingency_tables import mcnemar

from experiments.full_text.analyse_results import merge_predictions


def mcnemar_test(table):
    """Perform McNemar's test and return p-value."""
    result = mcnemar(table, exact=False, correction=True)
    return result.pvalue


classification_results_folder = "results/"
prompting_results_folder = "../../data/processed/prompting/"

csmed_ft = load_dataset(
    "../../ec2s/big_screening/datasets/csmed_ft", name="csmed_ft_all_source"
)


def get_predictions_classification(model_name, split):
    predictions_file = f"{classification_results_folder}/{model_name}_{split}.json"

    if not os.path.exists(predictions_file):
        print(f"Missing {predictions_file}")
        raise FileNotFoundError

    with open(predictions_file, "r") as f:
        predictions = json.load(f)

    y_pred = [np.argmax(pred) for pred in predictions]
    return y_pred


cls_models = [
    "yikuan8_Clinical-BigBird",
    "yikuan8_Clinical-Longformer",
    "allenai_longformer-base-4096",
    "google_bigbird-roberta-base",
]

gpt_models = [
    "predictions_gpt-3.5-turbo-16k.csv",
    "predictions_gpt-4-0314-50.csv",
    "predictions_gpt-3.5-turbo-0301-50.csv",
]


def contingency_table(pred_1, pred_2, gold):
    a = b = c = d = 0
    for i in range(len(pred_1)):
        if pred_1[i] == pred_2[i]:
            if pred_1[i] == gold[i]:
                a += 1
            else:
                d += 1
        else:
            if pred_1[i] == gold[i]:
                b += 1
            elif pred_2[i] == gold[i]:
                c += 1
    return np.array([[a, b], [c, d]])


if __name__ == "__main__":

    for split in ["sample", "test"]:
        dataset = csmed_ft[split]
        y_true = dataset["label"]

        combined_predictions = {}
        for model_name in cls_models:
            y_pred = get_predictions_classification(model_name, split)
            print(len(y_pred))
            print(f"Results for {model_name} on {split}, {sum(y_pred)} positive")
            combined_predictions[model_name] = y_pred

        if split == "sample":
            for model_name in gpt_models:
                predictions_file = f"{prompting_results_folder}{model_name}"
                predictions = pd.read_csv(predictions_file)
                predictions = merge_predictions(predictions)
                print(len(predictions))

                y_pred = predictions["prediction"]
                print(f"Results for {model_name} on {split}, {sum(y_pred)} positive")
                combined_predictions[model_name] = y_pred
        else:
            gpt_model = "test_predictions_gpt-3.5-turbo-16k.csv"
            predictions_file = f"{prompting_results_folder}{gpt_model}"
            predictions = pd.read_csv(predictions_file)
            predictions = merge_predictions(predictions)
            print(len(predictions))

            y_pred = predictions["prediction"]
            print(f"Results for {gpt_model} on {split}, {sum(y_pred)} positive")
            combined_predictions[gpt_model] = y_pred

        baseline_model = DummyClassifier(strategy="stratified", random_state=1)
        baseline_model.fit(np.zeros(len(y_true)), y_true)
        y_pred = baseline_model.predict(np.zeros(len(y_true)))
        print(f"Results for baseline on {split}, {sum(y_pred)} positive")
        combined_predictions["baseline"] = y_pred

        for model in combined_predictions:
            print(
                f"{model}: {sum(combined_predictions[model])} positive ({100 * sum(combined_predictions[model]) / len(combined_predictions[model])}%)"
            )

        results = {}
        tables = []
        for index_i, model_1 in enumerate(combined_predictions):
            for index_j, model_2 in enumerate(combined_predictions):
                if index_j > index_i:
                    print(f"Comparing {model_1} and {model_2}")
                    table = contingency_table(
                        combined_predictions[model_1],
                        combined_predictions[model_2],
                        y_true,
                    )
                    tables.append(table)
                    results[f"{model_1} vs {model_2}"] = table

        p_values = [mcnemar_test(table) for table in tables]

        alpha = 0.05
        bonferroni_corrected_alpha = alpha / len(tables)
        print(len(tables))

        print(f"Bonferroni corrected significance level: {bonferroni_corrected_alpha}")
        for i, p_value in enumerate(p_values):
            if p_value < bonferroni_corrected_alpha:
                print(
                    f"Significant difference between {list(results.keys())[i]} with p-value {p_value}"
                )
            else:
                print(
                    f"No significant difference between {list(results.keys())[i]} with p-value {p_value}"
                )
