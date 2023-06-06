import pandas as pd
import tiktoken
from datasets import load_dataset

enc = tiktoken.encoding_for_model("gpt-4")

livsb_ft = load_dataset(
    "../ec2s/big_screening/datasets/livsb_ft", name="livsb_ft_all_source"
)


def encode_livsb_ft(dataset):
    tokenised_dataset = []
    for example in dataset:
        review = (
            example["review_title"]
            + "\n"
            + example["review_abstract"]
            + "\n"
            + example["review_criteria"]
        ).strip()
        paper = (
            example["title"]
            + "\n\n"
            + example["abstract"]
            + "\n\n"
            + example["main_text"]
        ).strip()

        tokenised_dataset.append(
            {
                "review": len(enc.encode(review)),
                "paper": len(enc.encode(paper)),
                "total": len(enc.encode(review + "\n\n" + paper)),
            }
        )

    return tokenised_dataset


for split in ["train", "validation", "test", "sample"]:
    df = pd.DataFrame(encode_livsb_ft(livsb_ft[split]))
    print(split)
    print(df.describe())
    print()
