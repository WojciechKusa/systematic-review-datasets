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
import pickle
from typing import List, Tuple, Dict, Any

import datasets
import pandas as pd

from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.loader.bigbiohub import text_features

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
@article{@inproceedings{alharbi2019dataset,
  title={A dataset of systematic review updates},
  author={Alharbi, Amal and Stevenson, Mark},
  booktitle={Proceedings of the 42nd International ACM SIGIR Conference on Research and Development in Information Retrieval},
  pages={1257--1260},
  year={2019}
  }
"""

_DATASETNAME = "srupdates"
_DISPLAYNAME = "srupdates"

_DESCRIPTION = """\
Systematic Review Update Dataset - a dataset created to evaluate the retrieval performance in systematic reviews update.
"""

_HOMEPAGE = "https://github.com/Amal-Alharbi/Systematic_Reviews_Update"
_LICENSE = ""

_URLS = {
    "srupdates": "https://raw.githubusercontent.com/Amal-Alharbi/Systematic_Reviews_Update/master/update_dataset.pkl.zip",
}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION, Tasks.QUESTION_ANSWERING]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]
REVIEWS = [
    "CD000155",
    "CD000160",
    "CD000523",
    "CD001298",
    "CD001552",
    "CD002064",
    "CD002733",
    "CD004069",
    "CD004214",
    "CD004241",
    "CD004479",
    "CD005025",
    "CD005055",
    "CD005083",
    "CD005128",
    "CD005426",
    "CD005607",
    "CD006839",
    "CD006902",
    "CD007020",
    "CD007428",
    "CD008127",
    "CD008392",
    "CD010089",
    "CD010847",
]


class SrUpdatesDataset(datasets.GeneratorBasedBuilder):
    """Systematic Review Update Dataset - a dataset created to evaluate the retrieval performance in systematic reviews update.
    It consists of 25 Intervention systematic reviews."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    dataset_versions = ["all"]
    for dataset_version in dataset_versions:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"srupdates_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"srupdates {dataset_version} source schema",
                schema="source",
                subset_id=f"srupdates_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"srupdates_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"srupdates {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"srupdates_{dataset_version}",
            )
        )

    DEFAULT_CONFIG_NAME = "srupdates_all_source"

    def _info(self) -> datasets.DatasetInfo:

        if self.config.schema == "source":
            features = datasets.Features(
                {
                    "review_id": datasets.Value("string"),
                    "review_title": datasets.Value("string"),
                    "pmid": datasets.Value("string"),
                    "title": datasets.Value("string"),
                    "abstract": datasets.Value("string"),
                    "mesh_terms": [datasets.Value("string")],
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
        data_dir = dl_manager.download_and_extract(_URLS["srupdates"])

        reviews_pickle = pickle.load(open(f"{data_dir}/update_dataset.pkl", "rb"))

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                gen_kwargs={
                    "reviews_pickle": reviews_pickle,
                    "split": "train",
                },
            ),
        ]

    def _generate_examples(
        self,
        reviews_pickle: dict[str, Any],
        split: str,
    ) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        review = "_".join(self.config.subset_id.split("_")[1:])

        uid = 0

        _examples = []
        if review == "all":
            reviews = reviews_pickle.keys()
        else:
            reviews = [review]

        for r in reviews:
            for study in reviews_pickle[r]["search_results"]["original_records"]:
                if study[0] in reviews_pickle[r]["included"]["original_abs"]:
                    label = 1
                else:
                    label = 0
                _example = {
                    "review_id": r,
                    "review_title": reviews_pickle[r]["title"],
                    "pmid": study[0],
                    "title": study[1],
                    "abstract": study[2],
                    "mesh_terms": study[3],
                    "label": label,
                }
                _examples.append(_example)

        df = pd.DataFrame(_examples)

        for key, example in df.iterrows():

            title = example["title"]
            abstract = example["abstract"]
            label = example["label"]
            text = f"{title}\n\n{abstract}"

            uid += 1
            if self.config.schema == "source":
                data = {
                    "review_id": example["review_id"],
                    "review_title": example["review_title"],
                    "pmid": example["pmid"],
                    "title": title,
                    "abstract": abstract,
                    "mesh_terms": example["mesh_terms"],
                    "label": label,
                }
                yield str(uid), data

            elif self.config.schema == "bigbio_text":
                data = {
                    "id": str(uid),
                    "document_id": example["pmid"],
                    "text": text,
                    "labels": [label],
                }
                yield str(uid), data


if __name__ == "__main__":
    x = datasets.load_dataset(__file__, name="srupdates_all_source")
    print(type(x))
    print(x)
