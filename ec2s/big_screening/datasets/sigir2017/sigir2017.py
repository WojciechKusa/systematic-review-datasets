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
import os
import re
from typing import List, Tuple, Dict

import datasets
import pandas as pd

from ec2s.big_screening.datasets.sigir2017.prepare import prepare_dataset
from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.loader.bigbiohub import text_features

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
@inproceedings{scells2017collection,
	Author = {Scells, Harrisen and Zuccon, Guido and Koopman, Bevan and Deacon, Anthony and Geva, Shlomo and Azzopardi, Leif},
	Booktitle = {Proceedings of the 40th international ACM SIGIR conference on Research and development in Information Retrieval},
	Organization = {ACM},
	Title = {A Test Collection for Evaluating Retrieval of Studies for Inclusion in Systematic Reviews},
	Year = {2017}
}
"""

_DATASETNAME = "sigir2017"
_DISPLAYNAME = "sigir2017"

_DESCRIPTION = """\
Dataset containing a collection of queries for the paper "A Test Collection for Evaluating Retrieval of Studies for Inclusion in Systematic Reviews".
"""

_HOMEPAGE = "https://github.com/ielab/SIGIR2017-SysRev-Collection"
_LICENSE = ""

_URLS = {
    "sigir2017": "https://github.com/ielab/SIGIR2017-SysRev-Collection/archive/refs/heads/master.zip",
}
MAIN_REPO_FOLDER = "SIGIR2017-SysRev-Collection-master"

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION, Tasks.QUESTION_ANSWERING]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]


class Sigir2017Dataset(datasets.GeneratorBasedBuilder):
    """Dataset containing a collection of queries for the paper "A Test Collection for Evaluating Retrieval of Studies for Inclusion in Systematic Reviews"."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    dataset_versions = ["all"]
    for dataset_version in dataset_versions:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"sigir2017_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"sigir2017 {dataset_version} source schema",
                schema="source",
                subset_id=f"sigir2017_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"sigir2017_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"sigir2017 {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"sigir2017_{dataset_version}",
            )
        )

    DEFAULT_CONFIG_NAME = "sigir2017_all_source"

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
        data_dir = dl_manager.download_and_extract(_URLS["sigir2017"])
        pubmed_output_dir = "/".join(self.cache_dir.split("/")[:-3])

        input_folder = f"{data_dir}/{MAIN_REPO_FOLDER}"
        prepare_dataset(
            input_folder=input_folder,
            output_folder=pubmed_output_dir,
        )

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                gen_kwargs={
                    "qrels_data_dir": input_folder,
                    "docs_data_dir": pubmed_output_dir,
                    "split": "train",
                },
            ),
        ]

    def _generate_examples(
        self,
        qrels_data_dir: str,
        docs_data_dir: str,
        split: str,
    ) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        review = "_".join(self.config.subset_id.split("_")[1:])

        with open(f"{qrels_data_dir}/systematic_reviews.json") as f:
            reviews_mapping = json.load(f)
        cochrane_id_pattern = r"CD(\d+)"

        REVIEWS = [
            re.search(cochrane_id_pattern, r["url"]).group(0) for r in reviews_mapping
        ]
        uid = 0

        if review == "all":
            df = pd.DataFrame()
            for r in REVIEWS:
                review_df = pd.read_csv(os.path.join(docs_data_dir, f"{r}.csv"))
                review_df["Review"] = r
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
    x = datasets.load_dataset(__file__, name="sigir2017_all_source")
    print(type(x))
    print(x)
