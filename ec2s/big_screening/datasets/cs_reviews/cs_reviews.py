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
import json
from typing import List, Tuple, Dict

import bibtexparser
import datasets
import pandas as pd

from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.loader.bigbiohub import text_features

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
@inproceedings{hannousse2022semi,
  title={A Semi-automatic Document Screening System for Computer Science Systematic Reviews},
  author={Hannousse, Abdelhakim and Yahiouche, Salima},
  booktitle={Mediterranean Conference on Pattern Recognition and Artificial Intelligence},
  pages={201--215},
  year={2022},
  organization={Springer}
}
"""

_DATASETNAME = "cs_reviews"
_DISPLAYNAME = "cs_reviews"

_DESCRIPTION = """\
Seven systematic reviews from the computer science domain.
"""

_HOMEPAGE = "https://github.com/hannousse/Semantic-Scholar-Evaluation/"
_LICENSE = ""

_URLS = {
    "reviews_content": "https://github.com/hannousse/Semantic-Scholar-Evaluation/archive/refs/heads/master.zip",
    "reviews_metadata": "https://github.com/aliromagnoli/Erasmus-Stage/archive/refs/heads/main.zip",
}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]

REVIEWS = [
    "Alhammad-2018",
    "Ghasemi-2019",
    "Goulao-2016",
    "Guinea-2016",
    "Santos-2018",
    "Shahin-2017",
    "Yang-2016",
]


def prepare_cs_dataset(
    dataset_name: str,
    included_studies_repository: str,
    review_metadata_repository: str,
) -> pd.DataFrame:
    """
    :param dataset_name: name of the dataset
    :param included_studies_repository: path to the included studies repository
    :param review_metadata_repository: path to the review metadata repository
    :return:
    """

    included_studies_repository_path = f"{included_studies_repository}/Semantic-Scholar-Evaluation-master/Datasets for automatic screening of papers/{dataset_name}"

    included_file = f"{included_studies_repository_path}/{dataset_name}-Included.bib"
    excluded_file = f"{included_studies_repository_path}/{dataset_name}-Excluded.bib"
    review_metadata_file = f"{review_metadata_repository}/Erasmus-Stage-main/datasets/new_datasets/{dataset_name.lower()}.json"

    out_df = pd.DataFrame()
    for file, label in zip([included_file, excluded_file], [1, 0]):
        with open(file) as f:
            bib_database = bibtexparser.load(f)
            df = pd.DataFrame(bib_database.entries)
            df["label"] = label
            out_df = pd.concat([out_df, df], ignore_index=True)

    try:
        with open(review_metadata_file) as f:
            review_metadata = json.load(f)
    except json.decoder.JSONDecodeError:
        review_metadata = {}

    out_df["review_id"] = dataset_name
    return out_df


class CsReviewsDataset(datasets.GeneratorBasedBuilder):
    """Seven systematic reviews from the computer science domain."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []

    for dataset_version in REVIEWS:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"cs_reviews_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"cs_reviews {dataset_version} source schema",
                schema="source",
                subset_id=f"cs_reviews_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"cs_reviews_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"cs_reviews {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"cs_reviews_{dataset_version}",
            )
        )

    BUILDER_CONFIGS.append(
        BigBioConfig(
            name="cs_reviews_all_source",
            version=SOURCE_VERSION,
            description="cs_reviews all source schema",
            schema="source",
            subset_id="cs_reviews_all",
        )
    )

    DEFAULT_CONFIG_NAME = "cs_reviews_all_source"

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
        data_dir = dl_manager.download_and_extract(_URLS)

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                gen_kwargs={
                    "split": "train",
                    "data_dir": data_dir,
                },
            ),
        ]

    def _generate_examples(
        self,
        data_dir: dict[str, str],
        split: str,
    ) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        review = "_".join(self.config.subset_id.split("_")[2:])

        uid = 0

        if review == "all":
            reviews = REVIEWS
        else:
            reviews = [review]

        df = pd.DataFrame()
        for r in reviews:
            df = pd.concat(
                [
                    df,
                    prepare_cs_dataset(
                        dataset_name=r,
                        included_studies_repository=data_dir["reviews_content"],
                        review_metadata_repository=data_dir["reviews_metadata"],
                    ),
                ]
            )

        for key, example in df.iterrows():

            title = example["title"]
            abstract = example["abstract"]
            label = example["label"]
            review_id = example["review_id"]
            document_id = example["ID"]
            text = f"{title}\n\n{abstract}"

            uid += 1
            if self.config.schema == "source":
                data = {
                    "review_id": review_id,
                    "doc_id": document_id,
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


if __name__ == "__main__":
    x = datasets.load_dataset(__file__, name="cs_reviews_all_source")
    print(x)
