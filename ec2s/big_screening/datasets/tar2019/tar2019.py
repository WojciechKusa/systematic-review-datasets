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

from ec2s.big_screening.loader.bigbiohub import text_features
from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.datasets.tar2019.prepare import prepare_dataset


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


class Tar2019Dataset(datasets.GeneratorBasedBuilder):
    """Technologically Assisted Reviews in Empirical Medicine 2019."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    dataset_versions = ["all"]
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

    DEFAULT_CONFIG_NAME = "tar2019_source"

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

    def _generate_examples(
        self,
        qrels_data_dir: str,
        docs_data_dir: str,
        qrels_dict: dict[str, str],
        split: str,
    ) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        review = "_".join(self.config.subset_id.split("_")[1:])

        qrels_df = pd.DataFrame()
        for review_type, qrels_file in qrels_dict.items():
            qrels_path = os.path.join(qrels_data_dir, qrels_file)
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

        # breakpoint()
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
