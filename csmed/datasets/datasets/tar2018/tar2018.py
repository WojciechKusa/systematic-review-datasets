# coding=utf-8
# Copyright 2023 Wojciech Kusa
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os
from typing import List, Tuple, Dict

import datasets
import pandas as pd

from csmed.datasets.loader.bigbiohub import BigBioConfig
from csmed.datasets.loader.bigbiohub import Tasks
from csmed.datasets.loader.bigbiohub import text_features
from csmed.datasets.utils import (
    is_prepared,
    get_from_pubmed,
    save_checksum,
    mark_all_files_prepared,
)

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
@article{Kanoulas2018CLEFOverview,
	author = {Kanoulas, Evangelos and Li, Dan and Azzopardi, Leif and Spijker, Rene},
	issn = {1613-0073},
	journal = {CEUR Workshop Proceedings},
	keywords = {Cochrane, DTA, PubMed, TAR, active learning, benchmarking, diagnostic test accuracy, e-health, evaluation, high recall, information retrieval, relevance feedback, systematic reviews, technology assisted reviews, test collection, text classification},
	month = {7},
	title = {{CLEF 2018 technologically assisted reviews in empirical medicine overview}},
	url = {https://pureportal.strath.ac.uk/en/publications/clef-2018-technologically-assisted-reviews-in-empirical-medicine-},
	volume = {2125},
	year = {2018},
	bdsk-url-1 = {https://pureportal.strath.ac.uk/en/publications/clef-2018-technologically-assisted-reviews-in-empirical-medicine-}}
"""

_DATASETNAME = "tar2018"
_DISPLAYNAME = "tar2018"

_DESCRIPTION = """\
Technologically Assisted Reviews in Empirical Medicine 2018
"""

_HOMEPAGE = "https://github.com/CLEF-TAR/tar"
_LICENSE = "MIT license"

_URLS = {
    "tar": "https://github.com/WojciechKusa/tar/archive/refs/heads/master.zip",
}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]


def prepare_dataset(
    input_folder: str,
    output_folder: str,
    train_qrels: str,
    test_qrels: str,
) -> None:
    if is_prepared(output_folder):
        return

    qrels_df = pd.concat(
        [
            pd.read_csv(
                f"{input_folder}/{train_qrels}",
                sep="\s+",
                header=None,
                names=["review_id", "0", "PMID", "Label"],
            ),
            pd.read_csv(
                f"{input_folder}/{test_qrels}",
                sep="\s+",
                header=None,
                names=["review_id", "0", "PMID", "Label"],
            ),
        ]
    )

    print("PubMed data is being downloaded. This may take a while for the first time.")
    for review_id in qrels_df["review_id"].unique():
        review_df = qrels_df[qrels_df["review_id"] == review_id]
        print(f"{review_id=}, {len(review_df)=}")
        review_df = get_from_pubmed(review_df)
        review_df.to_csv(f"{output_folder}/{review_id}.csv", index=False)
        print(f"Prepared review size: {len(review_df)}")
        save_checksum(
            file=f"{output_folder}/{review_id}.csv", dataset_directory=output_folder
        )

    mark_all_files_prepared(output_folder)


class Tar2018Dataset(datasets.GeneratorBasedBuilder):
    """Technologically Assisted Reviews in Empirical Medicine 2018."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    dataset_versions = ["all"]
    for dataset_version in dataset_versions:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"tar2018_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"tar2018 {dataset_version} source schema",
                schema="source",
                subset_id=f"tar2018_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"tar2018_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"tar2018 {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"tar2018_{dataset_version}",
            )
        )

    DEFAULT_CONFIG_NAME = "tar2018_all_source"

    def _info(self) -> datasets.DatasetInfo:

        if self.config.schema == "source":
            features = datasets.Features(
                {
                    "review_name": datasets.Value("string"),
                    "pmid": datasets.Value("string"),
                    "title": datasets.Value("string"),
                    "abstract": datasets.Value("string"),
                    "label": datasets.ClassLabel(names=_CLASS_NAMES),
                }
            )
        elif self.config.schema == "bigbio_text":
            features = text_features
        else:
            raise ValueError(f"Unsupported schema {self.config.schema}")

        return datasets.DatasetInfo(
            description=_DESCRIPTION,
            features=features,
            homepage=_HOMEPAGE,
            license=_LICENSE,
            citation=_CITATION,
        )

    def _split_generators(self, dl_manager) -> List[datasets.SplitGenerator]:
        """Returns SplitGenerators."""
        data_dir = dl_manager.download_and_extract(_URLS["tar"])
        pubmed_output_dir = "/".join(self.cache_dir.split("/")[:-3])

        train_qrels = (
            "tar-master/2018-TAR/Task2/Training/qrels/full.train.abs.2018.qrels"
        )
        test_qrels = "tar-master/2018-TAR/Task2/Testing/qrels/full.test.abs.2018.qrels"

        prepare_dataset(
            input_folder=data_dir,
            output_folder=pubmed_output_dir,
            train_qrels=train_qrels,
            test_qrels=test_qrels,
        )

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                gen_kwargs={
                    "qrels_path": os.path.join(data_dir, train_qrels),
                    "split": "train",
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.TEST,
                gen_kwargs={
                    "qrels_path": os.path.join(data_dir, test_qrels),
                    "split": "test",
                },
            ),
        ]

    def _generate_examples(self, qrels_path, split: str) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        data_dir = "/".join(self.cache_dir.split("/")[:-3])
        review = "_".join(self.config.subset_id.split("_")[1:])
        qrels_df = pd.read_csv(
            qrels_path,
            sep="\s+",
            header=None,
            names=["review_id", "0", "PMID", "Label"],
        )
        REVIEWS = qrels_df["review_id"].unique().tolist()
        uid = 0

        if review == "all":
            df = pd.DataFrame()
            for r in REVIEWS:
                review_df = pd.read_csv(os.path.join(data_dir, f"{r}.csv"))
                review_df["Review"] = r
                review_df = review_df.drop(columns=["Label"])

                review_df = review_df.merge(
                    qrels_df,
                    left_on=["PMID", "Review"],
                    right_on=["PMID", "review_id"],
                    how="left",
                )

                df = pd.concat([df, review_df])
        else:
            df = pd.read_csv(os.path.join(data_dir, f"{review}.csv"))
            df["Review"] = review

        for key, example in df.iterrows():
            review_name = str(example["Review"])
            title = str(example["Title"])
            abstract = str(example["Abstract"])
            print(example["Label"], str(example["PMID"]))
            label = int(example["Label"])  # fixme soem labels are NaN
            pmid = str(example["PMID"])
            uid += 1
            text = f"{title}\n\n{abstract}"

            if self.config.schema == "source":
                data = {
                    "review_name": review_name,
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "label": label,
                }
                print(data)
                yield str(uid), data

            elif self.config.schema == "bigbio_text":
                data = {
                    "id": str(uid),
                    "document_id": pmid,
                    "text": text,
                    "labels": [label],
                }
                yield str(uid), data


if __name__ == "__main__":
    x = datasets.load_dataset(__file__, name="tar2018_all_source")
    print(type(x))
    print(x)
