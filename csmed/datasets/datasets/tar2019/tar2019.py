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
@article{Kanoulas2019CLEFOverview,
	author = {E. Kanoulas and Dan Li and Leif Azzopardi and Ren{\'e} Spijker},
	booktitle = {CLEF},
	title = {{CLEF 2019 Technology Assisted Reviews in Empirical Medicine Overview}},
	year = {2019}}
"""

_DATASETNAME = "tar2019"
_DISPLAYNAME = "tar2019"

_DESCRIPTION = """\
Technologically Assisted Reviews in Empirical Medicine 2019
"""

_HOMEPAGE = "https://github.com/CLEF-TAR/tar"
_LICENSE = "MIT license"

_URLS = {
    "tar": "https://github.com/WojciechKusa/tar/archive/refs/heads/master.zip",
}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION, Tasks.QUESTION_ANSWERING]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]

REVIEWS = [
    "CD012567",
    "CD008081",
    "CD010653",
    "CD002143",
    "CD011145",
    "CD010239",
    "CD011549",
    "CD001261",
    "CD009925",
    "CD010633",
    "CD012551",
    "CD009185",
    "CD008803",
    "CD011571",
    "CD010864",
    "CD011380",
    "CD008054",
    "CD008892",
    "CD012930",
    "CD011926",
    "CD009323",
    "CD011787",
    "CD009519",
    "CD010296",
    "CD011975",
    "CD011515",
    "CD011140",
    "CD012216",
    "CD012083",
    "CD010386",
    "CD008587",
    "CD005253",
    "CD009694",
    "CD012281",
    "CD012009",
    "CD012455",
    "CD009642",
    "CD012669",
    "CD011420",
    "CD010705",
    "CD012233",
    # "CD012668",
    "CD012661",
    "CD005139",
    "CD010276",
    "CD010542",
    "CD012521",
    "CD008874",
    "CD006715",
    "CD010019",
    "CD012120",
    "CD010213",
    "CD007868",
    "CD008686",
    "CD012347",
    "CD007867",
    "CD011768",
    "CD009551",
    "CD009135",
    "CD010438",
    "CD012165",
    "CD011977",
    "CD009372",
    "CD010409",
    "CD009944",
    "CD010657",
    "CD009579",
    "CD011912",
    "CD012599",
    "CD012164",
    "CD009786",
    "CD010038",
    "CD012768",
    "CD007427",
    "CD011126",
    "CD011602",
    "CD012080",
    "CD011053",
    "CD009263",
    "CD009069",
    "CD012223",
    "CD010526",
    "CD011436",
    "CD011431",
    "CD010173",
    "CD012019",
    "CD010778",
    "CD007394",
    "CD012010",
    "CD010339",
    "CD011686",
    "CD008760",
    "CD010753",
    "CD009044",
    "CD006468",
    "CD010355",
    "CD008759",
    "CD009647",
    "CD010558",
    "CD000996",
    "CD009020",
    "CD010502",
    "CD012069",
    "CD010680",
    "CD011558",
    "CD008018",
    "CD009591",
    "CD004414",
    "CD012179",
    "CD012342",
    "CD011134",
    "CD010023",
]


def prepare_dataset(
    input_folder: str, output_folder: str, dataset_splits: dict[str, dict[str, str]]
) -> None:
    if is_prepared(output_folder):
        print("PubMed data is already prepared.")
        return

    for dataset_split, review_types in dataset_splits.items():
        for review_type, qrels_file in review_types.items():
            qrels_df = pd.read_csv(
                f"{input_folder}/tar-master/2019-TAR/Task2/{dataset_split}/{review_type}/qrels/{qrels_file}",
                sep="\s+",
                header=None,
                names=["review_id", "0", "PMID", "Label"],
            )

            print(
                "PubMed data is being downloaded. This may take a while for the first time."
            )
            for review_id in qrels_df["review_id"].unique():
                review_df = qrels_df[qrels_df["review_id"] == review_id]
                print(f"{review_id=}, {len(review_df)=}")
                review_df = get_from_pubmed(review_df)
                review_df.to_csv(f"{output_folder}/{review_id}.csv", index=False)
                print(f"Prepared review size: {len(review_df)}")
                save_checksum(
                    file=f"{output_folder}/{review_id}.csv",
                    dataset_directory=output_folder,
                )

            mark_all_files_prepared(output_folder)


class Tar2019Dataset(datasets.GeneratorBasedBuilder):
    """Technologically Assisted Reviews in Empirical Medicine 2019."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    dataset_versions = REVIEWS
    for dataset_version in dataset_versions:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"tar2019_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"tar2019 {dataset_version} source schema",
                schema="source",
                subset_id=f"tar2019_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"tar2019_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"tar2019 {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"tar2019_{dataset_version}",
            )
        )

    BUILDER_CONFIGS.append(
        BigBioConfig(
            name="tar2019_all_source",
            version=SOURCE_VERSION,
            description="tar2019 all source schema",
            schema="source",
            subset_id="tar2019_all",
        )
    )
    DEFAULT_CONFIG_NAME = "tar2019_all_source"

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

        dataset_splits = {
            "Training": {
                "Intervention": "full.train.int.abs.2019.qrels",
                "DTA": "full.train.dta.abs.2019.qrels",
            },
            "Testing": {
                "Intervention": "full.test.intervention.abs.2019.qrels",
                "DTA": "full.test.dta.abs.2019.qrels",
                "Prognosis": "full.test.prognosis.abs.2019.qrels",
                "Qualitative": "full.test.qualitative.abs.2019.qrels",
            },
        }

        prepare_dataset(
            input_folder=data_dir,
            output_folder=pubmed_output_dir,
            dataset_splits=dataset_splits,
        )
        dataset_version = self.config.subset_id.split("_")[1]
        if dataset_version == "all":
            return [
                datasets.SplitGenerator(
                    name=datasets.Split.TRAIN,
                    gen_kwargs={
                        "qrels_data_dir": data_dir,
                        "docs_data_dir": pubmed_output_dir,
                        "qrels_dict": dataset_splits["Training"],
                        "split": "train",
                    },
                ),
                datasets.SplitGenerator(
                    name=datasets.Split.TEST,
                    gen_kwargs={
                        "qrels_data_dir": data_dir,
                        "docs_data_dir": pubmed_output_dir,
                        "qrels_dict": dataset_splits["Testing"],
                        "split": "test",
                    },
                ),
            ]
        else:
            return [
                datasets.SplitGenerator(
                    name=datasets.Split.TRAIN,
                    gen_kwargs={
                        "qrels_data_dir": data_dir,
                        "docs_data_dir": pubmed_output_dir,
                        "qrels_dict": dataset_splits["Training"],
                        "split": "train",
                    },
                )
            ]

    def _generate_examples(
        self,
        qrels_data_dir: str,
        docs_data_dir: str,
        qrels_dict: dict[str, str],
        split: str,
    ) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        if split == "train":
            dataset_split = "Training"
        elif split == "test":
            dataset_split = "Testing"
        else:
            raise ValueError(f"Unsupported split {split}")

        review = "_".join(self.config.subset_id.split("_")[1:])

        qrels_df = pd.DataFrame()
        for review_type, qrels_file in qrels_dict.items():
            qrels_path = os.path.join(
                qrels_data_dir,
                f"tar-master/2019-TAR/Task2/{dataset_split}/{review_type}/qrels/{qrels_file}",
            )
            qrels_df = pd.concat(
                [
                    qrels_df,
                    pd.read_csv(
                        qrels_path,
                        sep="\s+",
                        header=None,
                        names=["review_id", "0", "PMID", "Label"],
                    ),
                ]
            )

        REVIEWS = qrels_df["review_id"].unique().tolist()
        uid = 0

        if review == "all":
            df = pd.DataFrame()
            for r in REVIEWS:
                review_df = pd.read_csv(os.path.join(docs_data_dir, f"{r}.csv"))
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
            df = pd.read_csv(os.path.join(docs_data_dir, f"{review}.csv"))
            df["Review"] = review

        for key, example in df.iterrows():
            review_name = example["Review"]
            title = example["Title"]
            abstract = example["Abstract"]
            label = example["Label"]
            try:
                pmid = str(example["PMID"])
            except:
                pmid = "NA"  # some reviews don't have PMIDs
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
    x = datasets.load_dataset(__file__, name="tar2019_all_source")
    print(type(x))
    print(x)
