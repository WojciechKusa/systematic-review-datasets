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

from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.loader.bigbiohub import text_features

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
@article{wallace2010semi,
  title={Semi-automated screening of biomedical citations for systematic reviews},
  author={Wallace, Byron C and Trikalinos, Thomas A and Lau, Joseph and Brodley, Carla and Schmid, Christopher H},
  journal={BMC bioinformatics},
  volume={11},
  number={1},
  pages={1--11},
  year={2010},
  publisher={BioMed Central}
}
"""

_DATASETNAME = "sas2010"
_DISPLAYNAME = "sas2010"

_DESCRIPTION = """\
Three clinical systematic reviews
"""

_HOMEPAGE = "https://github.com/bwallace/citation-screening/"
_LICENSE = ""

_URLS = {
    "sas2010": "https://github.com/bwallace/citation-screening/archive/refs/heads/master.zip",
}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]

REVIEWS = [
    {"dataset_name": "proton_beam", "labels_filename": "proton_labeled_features"},
    {
        "dataset_name": "micro_nutrients",
        "labels_filename": "micro_labeled_features_only",
    },
    {"dataset_name": "copd", "labels_filename": "copd_labeled_terms_only"},
]


def prepare_clinical_dataset(
    dataset_name: str,
    labels_filename: str,
    repository_path: str,
) -> pd.DataFrame:
    """
    github.com/bwallace/citation-screening/tree/master/modeling/curious_snake/data

    :param dataset_name: name of the dataset
    :param labels_filename: filename that contains labels for that dataset
    :param repository_path: path to the root of bwallace/citation-screening/ repository
    :return:
    """
    labels_column: str = "Label"

    in_path = f"{repository_path}/modeling/curious_snake/data/{dataset_name}"
    labels_file = f"{in_path}/{labels_filename}"

    abstract_dir = f"{in_path}/Abstracts/"
    title_dir = f"{in_path}/Titles/"

    dataset = {}
    with open(labels_file) as fp:
        for content in fp.readlines():
            paper = {}
            doc_id, label = content.split(" ")[:2]
            if label == "-1":
                paper[labels_column] = 0
            else:
                paper[labels_column] = label

            dataset[doc_id] = paper

    for abstract_file in os.listdir(abstract_dir):
        if (
            not os.path.isdir(f"{abstract_dir}/{abstract_file}")
            and abstract_file in dataset
        ):
            doc_id = abstract_file.split("/")[-1]
            with open(f"{abstract_dir}/{abstract_file}") as fp:
                abstract = fp.readline()

            with open(f"{title_dir}/{abstract_file}") as fp:
                title = fp.readline()

            dataset[doc_id]["Title"] = f"{title.lower()}. {abstract.lower()}"
            dataset[doc_id]["Abstract"] = f"{title.lower()}. {abstract.lower()}"

    df = pd.DataFrame.from_dict(dataset).transpose()
    df["review_id"] = dataset_name

    return df


class Sas2010Dataset(datasets.GeneratorBasedBuilder):
    """Three clinical systematic reviews."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    dataset_versions = [review["dataset_name"] for review in REVIEWS]
    for dataset_version in dataset_versions:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"sas2010_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"sas2010 {dataset_version} source schema",
                schema="source",
                subset_id=f"sas2010_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"sas2010_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"sas2010 {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"sas2010_{dataset_version}",
            )
        )

    BUILDER_CONFIGS.append(
        BigBioConfig(
            name="sas2010_all_source",
            version=SOURCE_VERSION,
            description=f"sas2010 all source schema",
            schema="source",
            subset_id=f"sas2010_all",
        )
    )

    DEFAULT_CONFIG_NAME = "sas2010_all_source"

    def _info(self) -> datasets.DatasetInfo:

        if self.config.schema == "source":
            features = datasets.Features(
                {
                    "review_id": datasets.Value("string"),
                    "doc_id": datasets.Value("string"),
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
        data_dir = dl_manager.download_and_extract(_URLS["sas2010"])
        repository_main = f"{data_dir}/citation-screening-master"

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                gen_kwargs={
                    "split": "train",
                    "repository_path": repository_main,
                },
            ),
        ]

    def _generate_examples(
        self,
        repository_path: str,
        split: str,
    ) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        review = "_".join(self.config.subset_id.split("_")[1:])

        uid = 0

        if review == "all":
            reviews = [r["dataset_name"] for r in REVIEWS]
        else:
            reviews = [review]

        df = pd.DataFrame()
        for r in reviews:
            df = pd.concat(
                [
                    df,
                    prepare_clinical_dataset(
                        dataset_name=r,
                        labels_filename=[
                            x["labels_filename"]
                            for x in REVIEWS
                            if x["dataset_name"] == r
                        ][0],
                        repository_path=repository_path,
                    ),
                ]
            )

        for key, example in df.iterrows():

            title = example["Title"]
            abstract = example["Abstract"]
            label = example["Label"]
            review_id = example["review_id"]
            text = f"{title}\n\n{abstract}"

            uid += 1
            if self.config.schema == "source":
                data = {
                    "review_id": review_id,
                    "doc_id": str(uid),
                    "title": title,
                    "abstract": abstract,
                    "label": label,
                }
                yield str(uid), data

            elif self.config.schema == "bigbio_text":
                data = {
                    "id": str(uid),
                    "document_id": str(uid),
                    "text": text,
                    "labels": [label],
                }
                yield str(uid), data
