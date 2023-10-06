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

import datasets
import pandas as pd

from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.loader.bigbiohub import (
    text_features,
    entailment_features,
    pairs_features,
)

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
"""

_DATASETNAME = "csmed_ft"
_DISPLAYNAME = "csmed_ft"

_DESCRIPTION = """\
Systematic Review dataset focused on full-text screening. 
"""

_HOMEPAGE = "https://github.com/WojciechKusa/systematic-review-datasets"
_LICENSE = "CC BY-SA 4.0"

_URLS = {"csmed_ft": "../../../../data/CSMeD/CSMeD-FT.zip"}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION, Tasks.TEXTUAL_ENTAILMENT]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]


class CsmedFtDataset(datasets.GeneratorBasedBuilder):
    """Systematic Review dataset focused on full-text screening."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    dataset_versions = ["all"]
    for dataset_version in dataset_versions:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"csmed_ft_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"csmed_ft {dataset_version} source schema",
                schema="source",
                subset_id=f"csmed_ft_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"csmed_ft_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"csmed_ft {dataset_version} BigBio classification schema",
                schema="bigbio_text",
                subset_id=f"csmed_ft_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"csmed_ft_{dataset_version}_bigbio_te",
                version=BIGBIO_VERSION,
                description=f"csmed_ft {dataset_version} BigBio entailment schema",
                schema="bigbio_te",
                subset_id=f"csmed_ft_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"csmed_ft_{dataset_version}_bigbio_pairs",
                version=BIGBIO_VERSION,
                description=f"csmed_ft {dataset_version} BigBio text pairs classification schema",
                schema="bigbio_pairs",
                subset_id=f"csmed_ft_{dataset_version}",
            )
        )

    DEFAULT_CONFIG_NAME = "csmed_ft_all_source"

    def _info(self) -> datasets.DatasetInfo:

        if self.config.schema == "source":
            features = datasets.Features(
                {
                    "review_id": datasets.Value("string"),
                    "review_title": datasets.Value("string"),
                    "review_abstract": datasets.Value("string"),
                    "review_criteria": datasets.Value("string"),
                    "pmid": datasets.Value("string"),
                    "title": datasets.Value("string"),
                    "abstract": datasets.Value("string"),
                    "main_text": datasets.Value("string"),
                    "label": datasets.ClassLabel(names=_CLASS_NAMES),
                }
            )
        elif self.config.schema == "bigbio_text":
            features = text_features
        elif self.config.schema == "bigbio_te":
            features = entailment_features
        elif self.config.schema == "bigbio_pairs":
            features = pairs_features
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
        data_dir = dl_manager.download_and_extract(_URLS["csmed_ft"])

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                gen_kwargs={
                    "data_dir": data_dir,
                    "split": "train",
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.VALIDATION,
                gen_kwargs={
                    "data_dir": data_dir,
                    "split": "dev",
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.TEST,
                gen_kwargs={
                    "data_dir": data_dir,
                    "split": "test",
                },
            ),
            datasets.SplitGenerator(
                name="sample",
                gen_kwargs={
                    "data_dir": data_dir,
                    "split": "sample",
                },
            ),
        ]

    def _generate_examples(
        self,
        data_dir: str,
        split: str,
    ) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        documents_file = f"{data_dir}/CSMeD-FT/CSMeD-FT-{split}.csv"
        review_file = f"{data_dir}/CSMeD-FT/CSMeD-FT-{split}_reviews_metadata.json"

        df = pd.read_csv(documents_file)
        with open(review_file, "r") as f:
            reviews_metadata = json.load(f)

        uid = 0
        for key, example in df.iterrows():

            title = example["title"]
            abstract = example["abstract"]
            main_text = example["main_text"]
            label = example["decision"]
            if label == "included":
                label = 1
            elif label == "excluded":
                label = 0
            try:
                pmid = str(example["PubMed ID"])
            except ValueError:
                pmid = None

            review_criteria = reviews_metadata[example["review_id"]]["criteria_text"]
            review_abstract = reviews_metadata[example["review_id"]]["abstract"]
            review_title = reviews_metadata[example["review_id"]]["title"]

            text = f"{review_title}\n{review_abstract}\n{review_criteria}\n{title}\n\n{abstract}\n\n{main_text}"
            premise = f"{review_title}\n{review_abstract}\n{review_criteria}"
            hypothesis = f"{title}\n\n{abstract}\n\n{main_text}"

            uid += 1
            if self.config.schema == "source":
                data = {
                    "review_id": example["review_id"],
                    "review_title": review_title,
                    "review_abstract": review_abstract,
                    "review_criteria": review_criteria,
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "main_text": main_text,
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

            elif self.config.schema == "bigbio_te":
                data = {
                    "id": str(uid),
                    "premise": premise,
                    "hypothesis": hypothesis,
                    "label": [label],
                }
                yield str(uid), data
            elif self.config.schema == "bigbio_pairs":
                data = {
                    "id": str(uid),
                    "document_id": f"{example['review_id']}_{pmid}",  # combine review_id and pmid
                    "text_1": premise,
                    "text_2": hypothesis,
                    "label": [label],
                }
                yield str(uid), data
            else:
                raise ValueError(f"Unsupported schema {self.config.schema}")


if __name__ == "__main__":
    x = datasets.load_dataset(__file__, name="csmed_ft_all_source")
    print(type(x))
    print(x)
