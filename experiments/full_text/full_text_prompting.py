import os
from time import sleep

import openai
import pandas as pd
import tiktoken
from datasets import load_dataset
from langchain.text_splitter import TokenTextSplitter
from tqdm import tqdm

API_KEY_FILE = "../config/openai_key.txt"
with open(API_KEY_FILE, "r", encoding="utf-8") as f:
    api_key = f.read().strip()
openai.api_key = api_key


def take_only_n_words(text, n):
    return " ".join(text.split(" ")[:n])


def split_paper_into_chunks(paper, chunk_size):
    chunks = []
    for i in range(0, len(paper), chunk_size):
        chunks.append(paper[i : i + chunk_size])
    return chunks


def parse_response(response, example, uuid: int, paper_chunk_id: int) -> dict[str, str]:
    result = {
        "review_id": example["review_id"],
        "paper_id": example["pmid"],
        "response": response,
        "gold_label": example["label"],
        "uuid": uuid,
        "paper_chunk_id": paper_chunk_id,
    }
    answer = response["choices"][0].message.content
    if answer.lower() in ["included", "include", "yes"]:
        result["prediction"] = 1
    elif answer.lower() in ["excluded", "exclude", "no"]:
        result["prediction"] = 0
    else:
        result["prediction"] = -1

    return result


if __name__ == "__main__":
    outpath = "../data/processed/prompting/"

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    models = {"gpt-4-0314": 8000, "gpt-3.5-turbo-0301": 4000}

    livsb_ft = load_dataset(
        "../../ec2s/big_screening/datasets/livsb_ft", name="livsb_ft_all_source"
    )

    model_name = "gpt-3.5-turbo-0301"
    # model_name = "gpt-4-0314"
    context_size = models[model_name]
    review_splitter = TokenTextSplitter(
        chunk_size=context_size // 2, chunk_overlap=0, model_name=model_name
    )

    print(f"Using model {model_name}")
    enc = tiktoken.encoding_for_model(model_name)

    dataset = livsb_ft["sample"]

    predictions = []
    for index_i, example in tqdm(enumerate(dataset), total=len(dataset)):

        review = (
            "Systematic review:\n"
            + example["review_title"]
            + "\n"
            + example["review_criteria"]
        ).strip()
        paper = (
            "Candidate paper:\n"
            + example["title"]
            + "\n\n"
            + example["abstract"]
            + "\n\n"
            + example["main_text"]
        ).strip()

        len_review = len(enc.encode(review))
        paper_splitter = TokenTextSplitter(
            chunk_size=max(context_size // 2, context_size - len_review),
            chunk_overlap=0,
            model_name=model_name,
        )

        review = review_splitter.split_text(review)[0]
        paper_splits = paper_splitter.split_text(paper)
        print(f"Splitting paper into {len(paper_splits)} chunks")

        task_description = (
            "Does the following scientific paper fulfill all eligibility criteria and should it be included in the systematic review? "
            "Answer 'Included' or 'Excluded'."
        )
        for paper_chunk_id, paper in enumerate(paper_splits):
            prompt = (
                task_description
                + "\n\n"
                + review
                + "\n\n"
                + paper
                + "\n\n"
                + "Answer: "
            )
            print(len(enc.encode(prompt)))

            if model_name == "gpt-4-0314":
                while True:
                    try:
                        response = openai.ChatCompletion.create(
                            model=model_name,
                            messages=[
                                {
                                    "role": "user",
                                    "content": prompt,
                                },
                            ],
                        )
                        break
                    except Exception as e:
                        print(e)
                        print("Sleeping for 30 seconds")
                        sleep(30)
                        continue
            elif model_name == "gpt-3.5-turbo-0301":
                response = openai.ChatCompletion.create(
                    model=model_name,
                    messages=[
                        {
                            "role": "user",
                            "content": prompt,
                        },
                    ],
                )
            else:
                raise ValueError(f"Unknown model {model_name}")

            answer = parse_response(
                response, example=example, uuid=index_i, paper_chunk_id=paper_chunk_id
            )
            predictions.append(answer)

            df = pd.DataFrame(predictions)
            df.to_csv(f"{outpath}/predictions_{model_name}.csv", index=False)
