import json
import os

import evaluate
import numpy as np
import wandb
from datasets import load_dataset
from langchain.text_splitter import TokenTextSplitter
from transformers import AutoModelForSequenceClassification
from transformers import AutoTokenizer
from transformers import TrainingArguments, Trainer

wandb.init(project="full_text_classification")


def compute_metrics(eval_pred):
    metric1 = evaluate.load("precision")
    metric2 = evaluate.load("recall")
    metric3 = evaluate.load("f1")

    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)
    precision = metric1.compute(
        predictions=predictions, references=labels, average="macro"
    )["precision"]
    recall = metric2.compute(
        predictions=predictions, references=labels, average="macro"
    )["recall"]
    f1 = metric3.compute(predictions=predictions, references=labels)["f1"]
    f1_macro = metric3.compute(
        predictions=predictions, references=labels, average="macro"
    )["f1"]
    return {"precision": precision, "recall": recall, "f1": f1, "f1_macro": f1_macro}


def tokenize_function(examples):
    """Tokenize the text. We want the review to be the first sentence, then SEP token and then the paper."""
    reviews = [
        f"{review_title} {review_criteria}"
        for review_title, review_criteria in zip(
            examples["review_title"], examples["review_criteria"]
        )
    ]
    reviews = [review_splitter.split_text(review)[0] for review in reviews]
    papers = [
        f"{title} {abstract} {main_text}"
        for title, abstract, main_text in zip(
            examples["title"], examples["abstract"], examples["main_text"]
        )
    ]

    input_text = [
        review + tokenizer.sep_token + paper for review, paper in zip(reviews, papers)
    ]
    return tokenizer(input_text, padding="max_length", truncation=True)


if __name__ == "__main__":
    outpath = "../data/processed/prompting/"

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    models = {
        "allenai/longformer-base-4096": 4096,
        "yikuan8/Clinical-BigBird": 4096,
        "yikuan8/Clinical-Longformer": 4096,
    }

    dataset = load_dataset(
        "../../ec2s/big_screening/datasets/livsb_ft", name="livsb_ft_all_source"
    )
    for model_name, context_size in models.items():
        print(f"Using model {model_name}, context size {context_size}")

        review_splitter = TokenTextSplitter(
            chunk_size=context_size // 2,
            chunk_overlap=0,
            model_name="gpt-3.5-turbo-0301",
        )

        model = AutoModelForSequenceClassification.from_pretrained(
            model_name, num_labels=2
        )
        tokenizer = AutoTokenizer.from_pretrained(model_name)

        tokenized_datasets = dataset.map(tokenize_function, batched=True)

        training_args = TrainingArguments(
            output_dir="test_trainer",
            evaluation_strategy="epoch",
            per_device_train_batch_size=1,
            per_device_eval_batch_size=1,
            gradient_accumulation_steps=8,
            learning_rate=5e-5,
            num_train_epochs=20,
        )

        trainer = Trainer(
            model=model,
            args=training_args,
            train_dataset=tokenized_datasets["train"],
            eval_dataset=tokenized_datasets["validation"],
            compute_metrics=compute_metrics,
        )

        trainer.train()

        x = trainer.evaluate()
        print(x)

        for test_data in ["validation", "sample", "test"]:
            predictions = trainer.predict(
                tokenized_datasets[test_data],
            )
            # save predictions
            print(predictions)
            with open(f"{'_'.join(model_name.split('/'))}_{test_data}.json", "w") as fp:
                json.dump(predictions, fp)
