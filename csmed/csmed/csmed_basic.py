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

from typing import Union

import datasets

from csmed.loader.bigbiohub import Tasks
from csmed.loader.bigbiohub import text_features

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
"""

_DATASETNAME = "csmed_basic"
_DISPLAYNAME = "csmed_basic"

_DESCRIPTION = """\
30 systematic reviews from the biomedical and computer science domain.
Reviews contain list of papers with their titles, abstracts, and labels (included/excluded).
    """

_HOMEPAGE = "https://github.com/WojciechKusa/systematic-review-datasets"
_LICENSE = "CC BY-SA 4.0"

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]

CSMED_BASIC_REVIEWS = {
    "cohen": [
        "ACEInhibitors",
        "ADHD",
        "Antihistamines",
        "AtypicalAntipsychotics",
        "BetaBlockers",
        "CalciumChannelBlockers",
        "Estrogens",
        "NSAIDS",
        "Opiods",
        "OralHypoglycemics",
        "ProtonPumpInhibitors",
        "SkeletalMuscleRelaxants",
        "Statins",
        "Triptans",
        "UrinaryIncontinence",
    ],
    "cs_reviews": [
        "Alhammad-2018",
        "Ghasemi-2019",
        "Goulao-2016",
        "Guinea-2016",
        "Santos-2018",
        "Shahin-2017",
        "Yang-2016",
    ],
    "sas2010": [
        "proton_beam",
        "micro_nutrients",
        "copd",
    ],
    "swift": [
        "Neuropain",
        "Fluoride",
        "BPA",
        "Transgenerational",
        "PFOS-PFOA",
    ],
}


class CSMeDBasic:
    """30 systematic reviews from the biomedical and computer science domain.
    Reviews contain list of papers with their titles, abstracts, and labels (included/excluded).
    """

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    @staticmethod
    def load_dataset(
        base_path,
    ) -> dict[str, dict[str, dict[str, Union[str, datasets.Dataset]]]]:
        """Returns a dictionary of datasets.
        The keys are the names of the reviews.
        The values are the datasets.

        :param base_path: path to the directory where the CSMeDBasic file is stored"""
        csmed_basic_datasets = {}
        for dataset_name in CSMED_BASIC_REVIEWS:
            for review_name in CSMED_BASIC_REVIEWS[dataset_name]:
                _dataset = datasets.load_dataset(
                    path=f"{base_path}/../{dataset_name}/{dataset_name}.py",
                    name=f"{dataset_name}_{review_name}_bigbio_text",
                )
                csmed_basic_datasets[review_name] = {
                    "review_name": review_name,
                    "data": _dataset,
                }

        return {"TRAIN": csmed_basic_datasets}

    def _info(self) -> datasets.DatasetInfo:
        """Returns the dataset metadata."""

        features = text_features
        return datasets.DatasetInfo(
            description=_DESCRIPTION,
            features=features,
            homepage=_HOMEPAGE,
            license=_LICENSE,
            citation=_CITATION,
        )


if __name__ == "__main__":
    csmed_dataset = CSMeDBasic.load_dataset(base_path="")
    x = csmed_dataset.keys()
    print(x)
