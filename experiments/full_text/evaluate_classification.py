import json
import os

import numpy as np
import pandas as pd
from datasets import load_dataset

from experiments.full_text.analyse_results import get_evaluation

livsb_ft = load_dataset(
    "../../ec2s/big_screening/datasets/livsb_ft", name="livsb_ft_all_source"
)

for split in ["sample", "test"]:
    results = {}
    for model_name in [
        "yikuan8_Clinical-BigBird",
        "yikuan8_Clinical-Longformer",
        "allenai_longformer-base-4096",
        "google_bigbird-roberta-base",
    ]:

        predictions_file = f"results/{model_name}_{split}.json"

        if not os.path.exists(predictions_file):
            print(f"Missing {predictions_file}")
            continue

        with open(predictions_file, "r") as f:
            predictions = json.load(f)

        dataset = livsb_ft[split]
        y_pred = []
        for pred in predictions:
            y_pred.append(np.argmax(pred))
        y_true = dataset["label"]

        results[model_name] = get_evaluation(y_true, y_pred)

    print(f"Results for {split}")
    print(pd.DataFrame(results).T.round(3))
