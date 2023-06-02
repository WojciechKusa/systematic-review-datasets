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
import importlib

import datasets
import pandas as pd

from ec2s.big_screening.loader.bigbiohub import text_features
from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.datasets.tar2017.prepare import prepare_dataset


_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
@article{Kanoulas2018CLEFOverview,
	author = {Kanoulas, Evangelos and Li, Dan and Azzopardi, Leif and Spijker, Rene},
	issn = {1613-0073},
	journal = {CEUR Workshop Proceedings},
	keywords = {Cochrane, DTA, PubMed, TAR, active learning, benchmarking, diagnostic test accuracy, e-health, evaluation, high recall, information retrieval, relevance feedback, systematic reviews, technology assisted reviews, test collection, text classification},
	month = {7},
	title = {{CLEF 2018 technologically assisted reviews in empirical medicine overview}},
	url = {https://pureportal.strath.ac.uk/en/publications/clef-2018-technologically-assisted-reviews-in-empirical-medicine-},
	volume = {2125},
	year = {2018},
	bdsk-url-1 = {https://pureportal.strath.ac.uk/en/publications/clef-2018-technologically-assisted-reviews-in-empirical-medicine-}}
"""

_DATASETNAME = "tar2018"
_DISPLAYNAME = "tar2018"

_DESCRIPTION = """\
Technologically Assisted Reviews in Empirical Medicine 2018
"""

_HOMEPAGE = "https://github.com/CLEF-TAR/tar"
_LICENSE = "MIT license"

_URLS = {
    "tar": "https://github.com/WojciechKusa/tar/archive/refs/heads/master.zip",
}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"


_CLASS_NAMES = ["included", "excluded"]


class Tar2018Dataset(datasets.GeneratorBasedBuilder):
    """Technologically Assisted Reviews in Empirical Medicine 2018."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)


    BUILDER_CONFIGS = []
    dataset_versions = ["all"]
    for dataset_version in dataset_versions:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"tar2018_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"tar2018 {dataset_version} source schema",
                schema="source",
                subset_id=f"tar2018_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"tar2018_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"tar2018 {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"tar2018_{dataset_version}",
            )
        )

    DEFAULT_CONFIG_NAME = "tar2018_source"

    def _info(self) -> datasets.DatasetInfo:

        # Create the source schema; this schema will keep all keys/information/labels as close to the original dataset as possible.

        # You can arbitrarily nest lists and dictionaries.
        # For iterables, use lists over tuples or `datasets.Sequence`

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
        # TODO: This method is tasked with downloading/extracting the data and defining the splits depending on the configuration

        # If you need to access the "source" or "bigbio" config choice, that will be in self.config.name

        # LOCAL DATASETS: You do not need the dl_manager; you can ignore this argument. Make sure `gen_kwargs` in the return gets passed the right filepath

        # PUBLIC DATASETS: Assign your data-dir based on the dl_manager.

        # dl_manager is a datasets.download.DownloadManager that can be used to download and extract URLs; many examples use the download_and_extract method; see the DownloadManager docs here: https://huggingface.co/docs/datasets/package_reference/builder_classes.html#datasets.DownloadManager

        # dl_manager can accept any type of nested list/dict and will give back the same structure with the url replaced with the path to local files.
        dataset_name = self.config.name.split("_")[1]

        # TODO: KEEP if your dataset is PUBLIC; remove if not
        # urls = _URLS["tar"]
        data_dir = dl_manager.download_and_extract(_URLS["tar"])
        pubmed_output_dir = "/".join(self.cache_dir.split("/")[:-3])
        # data_dir = self.cache_dir
        # remove 'a41067f8dfcd8c869c22d97a584974a6235e420067f59c893f92b3509b80aba9' from '/Users/wojciechkusa/.cache/huggingface/datasets/bigbio_systematic_review/tar2018_cohen_source/1.0.0/a41067f8dfcd8c869c22d97a584974a6235e420067f59c893f92b3509b80aba9'
        # data_dir = "/".join(data_dir.split("/")[:-1])

        prepare_dataset(input_folder=data_dir, output_folder=pubmed_output_dir)

        # config_file = f"datasets/{dataset_name}/setup.json"
        # _preprocessing_module = (
        #     f"tar.big_screening.datasets.{dataset_name}.{dataset_name}"
        # )
        # prepare_dataset = getattr(
        #     importlib.import_module(_preprocessing_module), "prepare_dataset"
        # )

        # Not all datasets have predefined canonical train/val/test splits.
        # If your dataset has no predefined splits, use datasets.Split.TRAIN for all of the data.

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                # Whatever you put in gen_kwargs will be passed to _generate_examples
                gen_kwargs={
                    "qrels_path": os.path.join(data_dir, "tar-master/2018-TAR/Task2/Training/qrels/full.train.abs.2018.qrels"),
                    "split": "train",
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.TEST,
                gen_kwargs={
                    "qrels_path": os.path.join(data_dir, "tar-master/2018-TAR/Task2/Testing/qrels/full.test.abs.2018.qrels"),
                    "split": "test",
                },
            ),
            # datasets.SplitGenerator(
            #     name=datasets.Split.VALIDATION,
            #     gen_kwargs={
            #         "filepath": os.path.join(data_dir, "dev.jsonl"),
            #         "split": "dev",
            #     },
            # ),
        ]

    # method parameters are unpacked from `gen_kwargs` as given in `_split_generators`

    # TODO: change the args of this function to match the keys in `gen_kwargs`. You may add any necessary kwargs.

    def _generate_examples(self, qrels_path, split: str) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        data_dir = "/".join(self.cache_dir.split("/")[:-3])
        review = "_".join(self.config.subset_id.split("_")[1:])
        qrels_df = pd.read_csv(qrels_path, sep="\s+", header=None, names=["review_id", "0", "PMID", "Label"])
        REVIEWS = qrels_df["review_id"].unique().tolist()
        uid = 0
        # breakpoint()

        if review == "all":
            df = pd.DataFrame()
            for r in REVIEWS:
                review_df = pd.read_csv(os.path.join(data_dir, f"{r}.csv"))
                review_df["Review"] = r
                # remove old Label column
                review_df = review_df.drop(columns=["Label"])

                # add Label from qrels_df to review_df based on PMID and review_id
                review_df = review_df.merge(qrels_df, left_on=["PMID", "Review"], right_on=["PMID", "review_id"],
                                            how="left")

                df = pd.concat([df, review_df])
        else:
            df = pd.read_csv(os.path.join(data_dir, f"{review}.csv"))
            df["Review"] = review

        breakpoint()
        for key, example in df.iterrows():
            review_name = str(example["Review"])
            title = str(example["Title"])
            abstract = str(example["Abstract"])
            label = int(example["Label"])  # fixme soem labels are NaN
            pmid = str(example["PMID"])
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
                print(data)
                yield str(uid), data

            elif self.config.schema == "bigbio_text":
                data = {
                    "id": str(uid),
                    "document_id": pmid,
                    "text": text,
                    "labels": [label],
                }
                yield str(uid), data



        # TODO: This method handles input defined in _split_generators to yield (key, example) tuples from the dataset.

        # The `key` is for legacy reasons (tfds) and is not important in itself, but must be unique for each example.

        # # NOTE: For local datasets you will have access to self.config.data_dir and self.config.data_files
        #
        # if self.config.schema == "source":
        #     # TODO: yield (key, example) tuples in the original dataset schema
        #     for key, example in thing:
        #         yield key, example
        #
        # elif self.config.schema == "bigbio_text":
        #     # TODO: yield (key, example) tuples in the bigbio schema
        #     for key, example in thing:
        #         yield key, example
        #

# This template is based on the following template from the datasets package:
# https://github.com/huggingface/datasets/blob/master/templates/new_dataset_script.py


# This allows you to run your dataloader with `python tar.py` during development
# TODO: Remove this before making your PR
if __name__ == "__main__":
    x = datasets.load_dataset(__file__, name="tar2018_all_source")
    print(type(x))
    print(x)
