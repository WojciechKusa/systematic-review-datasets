import copy
import os

import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    fbeta_score,
)


def get_evaluation(y_true, y_pred):
    return {
        "acc": accuracy_score(y_true, y_pred),
        "prec": precision_score(y_true, y_pred, average="macro"),
        "recall": recall_score(y_true, y_pred, average="macro"),
        "f1": f1_score(y_true, y_pred, average="binary"),
        "f1-macro": f1_score(y_true, y_pred, average="macro"),
        "f1-micro": f1_score(y_true, y_pred, average="micro"),
        "f3": fbeta_score(y_true, y_pred, beta=3, average="macro"),
    }


def merge_predictions(df):
    """There are several prediction for the same record.
    Predictions can be 1, 0, or -1.
    take 1 if there are no 0.
    If there is at least one 0, take 0.
    If there are only -1, take 1.
    group by uuid"""
    new_df = copy.copy(df)
    new_df = new_df.groupby(["uuid"]).agg({"prediction": list}).reset_index()
    new_df["prediction"] = new_df["prediction"].apply(lambda x: [int(y) for y in x])
    new_df["prediction"] = new_df["prediction"].apply(lambda x: 0 if 0 in x else 1)

    new_df["gold_label"] = (
        df.groupby(["uuid"])
        .agg({"gold_label": list})
        .reset_index()["gold_label"]
        .apply(lambda x: x[0])
    )
    return new_df


if __name__ == "__main__":

    results_path = "../../data/processed/prompting/"

    model_predictions = [x for x in os.listdir(results_path) if x.endswith(".csv")]
    results = {}
    for model_name in model_predictions:
        df = pd.read_csv(f"{results_path}/{model_name}")
        df = merge_predictions(df)

        Y_true = df["gold_label"]
        Y_pred = df["prediction"]
        print(f"{model_name}, {len(df)} predictions")
        print(f"Number of included records: {sum(Y_true)}, {sum(Y_pred)}")
        print(
            f"there are {len(set(list(Y_true) + list(Y_pred)))} different labels in the dataset\n"
        )
        results["_".join(model_name.split("_")[1:])] = get_evaluation(Y_true, Y_pred)

    print(pd.DataFrame(results).T.round(3))
