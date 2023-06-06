import os

import evaluate
import numpy as np
from datasets import load_dataset
from langchain.text_splitter import TokenTextSplitter
from transformers import AutoModelForSequenceClassification
from transformers import AutoTokenizer
from transformers import TrainingArguments, Trainer


def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)
    return metric.compute(predictions=predictions, references=labels)


def tokenize_function(examples):
    return tokenizer(examples["main_text"], padding="max_length", truncation=True)


if __name__ == "__main__":
    outpath = "../data/processed/prompting/"

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    models = [
        "yikuan8/Clinical-BigBird",
        "yikuan8/Clinical-Longformer",
    ]
    context_size = 4096

    dataset = load_dataset(
        "../../ec2s/big_screening/datasets/livsb_ft", name="livsb_ft_all_source"
    )
    review_splitter = TokenTextSplitter(
        chunk_size=context_size // 2, chunk_overlap=0, model_name="gpt-3.5-turbo-0301"
    )
    model_name = models[1]
    print(f"Using model {model_name}")

    model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)
    tokenizer = AutoTokenizer.from_pretrained(model_name)

    tokenized_datasets = dataset.map(tokenize_function, batched=True)

    training_args = TrainingArguments(
        output_dir="test_trainer", evaluation_strategy="epoch"
    )

    metric = evaluate.load("precision")

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized_datasets["sample"],
        eval_dataset=tokenized_datasets["validation"],
        compute_metrics=compute_metrics,
    )

    trainer.train()

    trainer.evaluate()
