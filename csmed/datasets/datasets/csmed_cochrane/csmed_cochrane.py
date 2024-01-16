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
from typing import Union

import datasets

from csmed.datasets.loader.bigbiohub import Tasks
from csmed.datasets.loader.bigbiohub import text_features
from . import EVAL_REVIEWS, TRAIN_REVIEWS

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
"""

_DATASETNAME = "csmed_cochrane"
_DISPLAYNAME = "csmed_cochrane"

_DESCRIPTION = """\
300 medical systematic literature reviews produced by Cochrane.
Reviews contain list of papers with their titles, abstracts, and labels (included/excluded).
Reviews also contain its metadata (review abstract and eligibility criteria).
    """

_HOMEPAGE = "https://github.com/WojciechKusa/systematic-review-datasets"
_LICENSE = "CC BY-SA 4.0"

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]

CSMED_COCHRANE_REVIEWS = {
    "sigir2017": [
        "CD009784",
        "CD010206",
        "CD002764",
        "CD011146",
        "CD010264",
        "CD011721",
        "CD009782",
        "CD009771",
        "CD006962",
        "CD007470",
        "CD009946",
        "CD006569",
        "CD006900",
        "CD010910",
        "CD007953",
        "CD003504",
        "CD003160",
        "CD011209",
        "CD007505",
        "CD009400",
        "CD006014",
        "CD003559",
        "CD001800",
        "CD010375",
        "CD010381",
        "CD004920",
        "CD011066",
        "CD007103",
        "CD003557",
        "CD009436",
        "CD004944",
        "CD004178",
        "CD007736",
        "CD005059",
        "CD010927",
        "CD011837",
        "CD003353",
        "CD009669",
        "CD001457",
        "CD010387",
        "CD009634",
        "CD009633",
        "CD011241",
        "CD005019",
        "CD002832",
        "CD004903",
        "CD003914",
        "CD002062",
        "CD002250",
        "CD007379",
        "CD011472",
        "CD009210",
        "CD010139",
        "CD010333",
        "CD009887",
        "CD006832",
        "CD008963",
        "CD005614",
        "CD010767",
        "CD006748",
        "CD007497",
        "CD008873",
        "CD006770",
        "CD010411",
        "CD008482",
        "CD001957",
        "CD006172",
        "CD010225",
        "CD010473",
        "CD002783",
        "CD007037",
        "CD008881",
        "CD006127",
        "CD011732",
        "CD002117",
        "CD009593",
        "CD000284",
        "CD010685",
        "CD003266",
        "CD008265",
        "CD010632",
        "CD010269",
        "CD004250",
        "CD009780",
        "CD005340",
        "CD008807",
        "CD004054",
        "CD009125",
        "CD011724",
        "CD010266",
        "CD007028",
        "CD000031",
        "CD008435",
        "CD008056",
        "CD000009",
        "CD005346",
        "CD010693",
        "CD006995",
        "CD010235",
        "CD007271",
        "CD000313",
        "CD009403",
        "CD008345",
        "CD003137",
        "CD008170",
        "CD002283",
        "CD005260",
        "CD006612",
        "CD002073",
        "CD008185",
        "CD006678",
        "CD007356",
        "CD001431",
        "CD004520",
        "CD009831",
        "CD004376",
        "CD002840",
        "CD000384",
        "CD010912",
        "CD007719",
        "CD004396",
        "CD006638",
        "CD007525",
        "CD010901",
        "CD011447",
        "CD002898",
        "CD012004",
        "CD003344",
        "CD007115",
        "CD010709",
        "CD008366",
        "CD007716",
        "CD002069",
        "CD003987",
        "CD009685",
        "CD009823",
        "CD001843",
        "CD002201",
        "CD008969",
        "CD011047",
        "CD007122",
        "CD011281",
        "CD010533",
        "CD010534",
        "CD007745",
        "CD001150",
        "CD010479",
        "CD008812",
        "CD003855",
        "CD005397",
        "CD010226",
        "CD003067",
        "CD008020",
        "CD002115",
        "CD000014",
        "CD007699",
        "CD006546",
        "CD008426",
        "CD006745",
        "CD006342",
        "CD008213",
        "CD004288",
        "CD004875",
        "CD010447",
        "CD010845",
        "CD009530",
        "CD009537",
        "CD002143",
        "CD010425",
        "CD011145",
    ],
    "tar2019": [
        "CD012567",
        "CD008081",
        "CD010653",
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
    ],
    "sr_updates": [
        "CD004069",
        "CD008127",
        "CD001298",
        "CD010847",
        "CD004241",
        "CD010089",
        "CD007020",
        "CD006902",
        "CD005128",
        "CD007428",
        "CD002733",
        "CD005055",
        "CD008392",
        "CD006839",
        "CD005083",
        "CD004479",
    ],
}


class CSMeDCochrane:
    """296 systematic reviews from the biomedical and computer science domain.
    Reviews contain list of papers with their titles, abstracts, and labels (included/excluded).
    """

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    def load_dataset(
        self, base_path
    ) -> dict[str, dict[str, dict[str, Union[str, datasets.Dataset]]]]:
        """Returns a dictionary of datasets.
        The keys are the names of the reviews.
        The values are the datasets.

        :param base_path: path to the directory where the CSMeDCochrane file is stored
        """
        path = f"{base_path}/../../../../data/external/"
        csmed_cochrane_datasets = {"TRAIN": {}, "EVAL": {}}
        for dataset_name in CSMED_COCHRANE_REVIEWS:
            with open(f"{path}/{dataset_name}/data_index.json") as f:
                collection_details = json.load(f)

            for review_name in CSMED_COCHRANE_REVIEWS[dataset_name]:
                details_file = (
                    f"{path}/{dataset_name}/{review_name}/review_details.json"
                )
                with open(details_file) as f:
                    dataset_details = json.load(f)

                _dataset = datasets.load_dataset(
                    path=f"{base_path}/../{dataset_name}/{dataset_name}.py",
                    name=f"{dataset_name}_{review_name}_bigbio_text",
                )
                if review_name in TRAIN_REVIEWS:
                    csmed_cochrane_datasets["TRAIN"][review_name] = {
                        "review_name": review_name,
                        "dataset_details": dataset_details,
                        "data": _dataset,
                    }
                elif review_name in EVAL_REVIEWS:
                    csmed_cochrane_datasets["EVAL"][review_name] = {
                        "review_name": review_name,
                        "dataset_details": dataset_details,
                        "data": _dataset,
                    }
                else:
                    raise ValueError(f"Unknown review name: {review_name}")

        return csmed_cochrane_datasets

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
    csm = CSMeDCochrane()
    csmed_dataset = csm.load_dataset(base_path=".")

    x = csmed_dataset.keys()

    print(x)
