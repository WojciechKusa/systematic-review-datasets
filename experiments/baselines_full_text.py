import json
from tqdm import tqdm
import os
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, fbeta_score
from sklearn.dummy import DummyClassifier

data_path = "../data/LivSB/LivSB-FT/"
datasets = [x for x in os.listdir(data_path) if x.endswith('.csv')]


classifiers = {
    "most_frequent": DummyClassifier(strategy="most_frequent"),
    "stratified": DummyClassifier(strategy="stratified"),
    "uniform": DummyClassifier(strategy="uniform"),
    "constant_1": DummyClassifier(strategy="constant", constant=1),
}


for dataset in datasets:
    results = {}
    df = pd.read_csv(f"{data_path}/{dataset}")
    print(f"Dataset: {dataset}")
    print(f"Number of documents: {len(df)}")
    Y = df['decision']
    Y = Y.replace({'included': 1, 'excluded': 0})
    X = df['main_text']
    for name, clf in classifiers.items():
        clf.fit(X, Y)
        predictions = clf.predict(X)
        results[name] = {
            "accuracy": accuracy_score(Y, predictions),
            "precision": precision_score(Y, predictions),
            "recall": recall_score(Y, predictions),
            "f1": f1_score(Y, predictions),
            "f3": fbeta_score(Y, predictions, beta=3)
        }

    print(pd.DataFrame(results).T)