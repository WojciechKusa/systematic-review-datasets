import os

import evaluate
import numpy as np
from datasets import load_dataset
from langchain.text_splitter import TokenTextSplitter
from transformers import AutoModelForSequenceClassification
from transformers import AutoTokenizer
from transformers import TrainingArguments, Trainer

context_size = 4096
review_splitter = TokenTextSplitter(
    chunk_size=context_size // 2, chunk_overlap=0, model_name="gpt-3.5-turbo-0301"
)


def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)
    return metric.compute(predictions=predictions, references=labels)


def tokenize_function(examples):
    """ Tokenize the text. We want the review to be the first sentence, then SEP token and then the paper."""
    reviews = [f"{review_title} {review_criteria}" for review_title, review_criteria in zip(examples["review_title"], examples["review_criteria"])]
    reviews = [review_splitter.split_text(review)[0] for review in reviews]
    papers = [f"{title} {abstract} {main_text}" for title, abstract, main_text in zip(examples["title"], examples["abstract"], examples["main_text"])]

    input_text = [review + tokenizer.sep_token + paper for review, paper in zip(reviews, papers)]
    return tokenizer(input_text, padding="max_length", truncation=True)


if __name__ == "__main__":
    outpath = "../data/processed/prompting/"

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    models = [
        "yikuan8/Clinical-BigBird",
        "yikuan8/Clinical-Longformer",
    ]

    dataset = load_dataset(
        "../../ec2s/big_screening/datasets/livsb_ft", name="livsb_ft_all_source"
    )
    model_name = models[1]
    print(f"Using model {model_name}")

    model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)
    tokenizer = AutoTokenizer.from_pretrained(model_name)

    tokenized_datasets = dataset.map(tokenize_function, batched=True)
    metric = evaluate.load("precision")

    training_args = TrainingArguments(
        output_dir="test_trainer",
        evaluation_strategy="epoch",
        per_device_train_batch_size=2,
        per_device_eval_batch_size=2,
        gradient_accumulation_steps=8,
        learning_rate=3e-5,
        num_train_epochs=2,
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized_datasets["validation"],
        eval_dataset=tokenized_datasets["test"],
        compute_metrics=compute_metrics,
    )

    trainer.train()

    trainer.evaluate()

    # evaluate     tokenized_datasets['test']

