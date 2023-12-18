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
import copy
import json
import os
import re
from typing import List, Tuple, Dict

import datasets
import pandas as pd

from csmed.datasets.loader.bigbiohub import BigBioConfig
from csmed.datasets.loader.bigbiohub import Tasks
from csmed.datasets.loader.bigbiohub import text_features
from csmed.datasets.utils import (
    is_prepared,
    temporal_search_pubmed,
    get_from_pubmed,
    save_checksum,
    mark_all_files_prepared,
)

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
_LICENSE = "Unknown"

_URLS = {
    "sigir2017": "https://github.com/ielab/SIGIR2017-SysRev-Collection/archive/refs/heads/master.zip",
}
MAIN_REPO_FOLDER = "SIGIR2017-SysRev-Collection-master"

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION, Tasks.QUESTION_ANSWERING]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]

REVIEWS = [
    "CD010653",
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
    "CD002143",
    "CD010685",
    "CD003266",
    "CD011145",
    "CD008265",
    "CD009185",
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
    "CD010633",
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
    "CD001261",
    "CD006612",
    "CD010542",
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
    "CD008081",
    "CD010425",
]


class Sigir2017Dataset(datasets.GeneratorBasedBuilder):
    """Dataset containing a collection of queries for the paper "A Test Collection for Evaluating Retrieval of
    Studies for Inclusion in Systematic Reviews"."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    dataset_versions = REVIEWS
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

    BUILDER_CONFIGS.append(
        BigBioConfig(
            name="sigir2017_all_source",
            version=SOURCE_VERSION,
            description="sigir2017 all source schema",
            schema="source",
            subset_id=f"sigir2017_all",
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
cochrane_id_pattern = r"CD(\d+)"


def prepare_dataset(input_folder: str, output_folder: str) -> None:
    if is_prepared(output_folder):
        print("PubMed data is already prepared.")
        return

    qrels_df = pd.read_csv(
        f"{input_folder}/sigir2017.qrels",
        sep="\s+",
        header=None,
        names=["review_id", "0", "PMID", "Label"],
    )
    with open(f"{input_folder}/queries_unannotated.json") as f:
        queries = json.load(f)
    with open(f"{input_folder}/systematic_reviews.json") as f:
        reviews_mapping = json.load(f)

    print("PubMed data is being downloaded. This may take a while for the first time.")
    for review_id in qrels_df["review_id"].unique():  # review_ids are just numbers
        review_url = [r["url"] for r in reviews_mapping if r["id"] == review_id][0]
        match = re.search(cochrane_id_pattern, review_url)
        if match:
            cochrane_id = match.group(0)
        else:
            raise ValueError("Review ID not found")

        query = [q for q in queries if q["document_id"] == review_id]
        if len(query) > 0:
            n_found, id_list = temporal_search_pubmed(
                query[0]["query"], query[0]["start_date"], query[0]["end_date"]
            )
        else:
            n_found = 0
            id_list = []

        review_df = copy.copy(qrels_df[qrels_df["review_id"] == review_id])
        print(f"{cochrane_id=}, {len(review_df)=}, {len(id_list)=}")
        if len(id_list) > 0:
            review_df = pd.concat(
                [
                    review_df,
                    pd.DataFrame(
                        {
                            "review_id": review_id,
                            "0": 0,
                            "PMID": id_list,
                            "Label": 0,
                        }
                    ),
                ],
                ignore_index=True,
            )
        review_df["PMID"] = review_df["PMID"].astype(int)
        review_df = review_df.drop_duplicates(subset=["PMID"])

        review_df = get_from_pubmed(review_df)

        review_df["review_id"] = cochrane_id
        review_df["Label"] = review_df["Label"].apply(
            lambda x: 0 if x in [0, 1, 3] else 1
        )

        review_df.to_csv(f"{output_folder}/{cochrane_id}.csv", index=False)
        print(f"Prepared review size: {len(review_df)}")
        save_checksum(
            file=f"{output_folder}/{cochrane_id}.csv",
            dataset_directory=output_folder,
        )

    mark_all_files_prepared(output_folder)
