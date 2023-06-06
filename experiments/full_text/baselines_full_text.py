import os

import pandas as pd
from sklearn.dummy import DummyClassifier
from tqdm import tqdm

from experiments.full_text.analyse_results import get_evaluation

data_path = "../data/LivSB/LivSB-FT/"
datasets = [x for x in os.listdir(data_path) if x.endswith('.csv')]

random_states = list(range(1, 101))

classifiers = {
    "most_frequent": DummyClassifier(strategy="most_frequent"),
    "stratified": DummyClassifier(strategy="stratified"),
    "uniform": DummyClassifier(strategy="uniform"),
    "constant_0": DummyClassifier(strategy="constant", constant=0),
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
        if name in ["stratified", "uniform"]:
            averaged_results = {}
            for random_state in tqdm(random_states):
                clf = DummyClassifier(strategy=name, random_state=random_state)
                clf.fit(X, Y)
                predictions = clf.predict(X)
                averaged_results[random_state] = get_evaluation(Y, predictions)
            averaged_results = pd.DataFrame(averaged_results).T
            averaged_results = averaged_results.mean().to_dict()
            results[name] = averaged_results
        else:
            clf.fit(X, Y)
            predictions = clf.predict(X)

            results[name] = get_evaluation(Y, predictions)

    print(pd.DataFrame(results).T.round(3))
