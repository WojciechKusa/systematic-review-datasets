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

from ec2s.big_screening.datasets.cohen.prepare import prepare_dataset, REVIEWS
from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.loader.bigbiohub import text_features

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
@article{Cohen2006,
	author = {Cohen, A. M. and Hersh, W. R. and Peterson, K. and Yen, Po Yin},
	doi = {10.1197/jamia.M1929},
	issn = {10675027},
	journal = {Journal of the American Medical Informatics Association},
	month = {3},
	number = {2},
	pages = {206--219},
	pmid = {16357352},
	publisher = {Oxford University Press},
	title = {{Reducing workload in systematic review preparation using automated citation classification}},
	url = {/pmc/articles/PMC1447545/ /pmc/articles/PMC1447545/?report=abstract https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1447545/},
	volume = {13},
	year = {2006},
	bdsk-url-1 = {/pmc/articles/PMC1447545/%20/pmc/articles/PMC1447545/?report=abstract%20https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1447545/},
	bdsk-url-2 = {https://doi.org/10.1197/jamia.M1929}}
"""

_DATASETNAME = "Cohen"
_DISPLAYNAME = "Cohen"

_DESCRIPTION = """\
Systematic Drug Class Review Gold Standard Data
PubMed Identifiers Annotated by Inclusion in Systematic Review
Here is the data used in our research on automated classification of document citations for systematic review of drug classes.
"""

_HOMEPAGE = "https://dmice.ohsu.edu/cohenaa/systematic-drug-class-review-data.html"
_LICENSE = ""

_URLS = {
    _DATASETNAME: {
        "cohen": "https://dmice.ohsu.edu/cohenaa/epc-ir-data/epc-ir.clean.tsv",
    }
}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]


class CohenDataset(datasets.GeneratorBasedBuilder):
    """Systematic Drug Class Review Gold Standard Data."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    reviews = REVIEWS
    for dataset_version in reviews:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"cohen_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"cohen {dataset_version} source schema",
                schema="source",
                subset_id=f"cohen_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"cohen_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"cohen {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"cohen_{dataset_version}",
            )
        )

    # Add an "all" config that combines all the reviews -- only for source schema
    BUILDER_CONFIGS.append(
        BigBioConfig(
            name="cohen_all_source",
            version=SOURCE_VERSION,
            description="cohen all source schema",
            schema="source",
            subset_id="cohen_all",
        )
    )

    DEFAULT_CONFIG_NAME = "cohen_all_source"

    def _info(self) -> datasets.DatasetInfo:
        """Returns the dataset metadata."""

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

        data_dir = "/".join(self.cache_dir.split("/")[:-3])
        prepare_dataset(output_folder=data_dir)

        # Not all datasets have predefined canonical train/val/test splits.
        # If your dataset has no predefined splits, use datasets.Split.TRAIN for all of the data.

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                # Whatever you put in gen_kwargs will be passed to _generate_examples
                gen_kwargs={
                    "filepath": os.path.join(data_dir, "train.jsonl"),
                    "split": "train",
                },
            ),
        ]

    def _generate_examples(self, filepath, split: str) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        data_dir = "/".join(self.cache_dir.split("/")[:-3])
        review = "_".join(self.config.subset_id.split("_")[1:])

        uid = 0

        if review == "all":
            df = pd.DataFrame()
            for r in REVIEWS:
                review_df = pd.read_csv(os.path.join(data_dir, f"{r}.tsv"), sep="\t")
                review_df["Review"] = r
                df = pd.concat([df, review_df])
        else:
            df = pd.read_csv(os.path.join(data_dir, f"{review}.tsv"), sep="\t")
            df["Review"] = review

        for key, example in df.iterrows():
            review_name = example["Review"]
            title = example["Title"]
            abstract = example["Abstract"]
            label = example["Label"]
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
    x = datasets.load_dataset(__file__, name="cohen_ADHD_source")
    print(x)

    y = datasets.load_dataset(__file__, name="cohen_ADHD_bigbio_text")
    print(y)

    z = datasets.load_dataset(__file__, name="cohen_all_source")
    print(z)
