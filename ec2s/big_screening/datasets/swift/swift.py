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

from ec2s.big_screening.datasets.swift.prepare import prepare_dataset, REVIEWS
from ec2s.big_screening.loader.bigbiohub import BigBioConfig
from ec2s.big_screening.loader.bigbiohub import Tasks
from ec2s.big_screening.loader.bigbiohub import text_features

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
@article{Howard2016,
	author = {Howard, Brian E. and Phillips, Jason and Miller, Kyle and Tandon, Arpit and Mav, Deepak and Shah, Mihir R. and Holmgren, Stephanie and Pelch, Katherine E. and Walker, Vickie and Rooney, Andrew A. and Macleod, Malcolm and Shah, Ruchir R. and Thayer, Kristina},
	doi = {10.1186/s13643-016-0263-z},
	issn = {20464053},
	journal = {Systematic Reviews},
	keywords = {Literature prioritization, SWIFT-Review, Scoping reports, Software, Systematic review},
	month = {5},
	number = {1},
	pages = {1--16},
	pmid = {27216467},
	publisher = {BioMed Central Ltd.},
	title = {{SWIFT-Review: A text-mining workbench for systematic review}},
	url = {https://link.springer.com/articles/10.1186/s13643-016-0263-z https://link.springer.com/article/10.1186/s13643-016-0263-z},
	volume = {5},
	year = {2016},
	bdsk-url-1 = {https://link.springer.com/articles/10.1186/s13643-016-0263-z%20https://link.springer.com/article/10.1186/s13643-016-0263-z},
	bdsk-url-2 = {https://doi.org/10.1186/s13643-016-0263-z}}
"""

_DATASETNAME = "swift"
_DISPLAYNAME = "swift"

_DESCRIPTION = """\
Four datasets (Additional file 1) were generated by the National Toxicology Program (NTP) 
    Office of Health Assessment and Translation (OHAT), one dataset (Additional file 2) was provided by 
    the Edinburgh CAMARADES group (www.camarades.info).
    """

_HOMEPAGE = "https://systematicreviewsjournal.biomedcentral.com/articles/10.1186/s13643-016-0263-z#Sec30"
_LICENSE = "CC0 1.0 Universal (CC0 1.0) Public Domain Dedication"

_URLS = {
    _DATASETNAME: {
        "ohat": "https://static-content.springer.com/esm/art%3A10.1186%2Fs13643-016-0263-z/MediaObjects/13643_2016_263_MOESM1_ESM.xlsx",
        "camrades": "https://static-content.springer.com/esm/art%3A10.1186%2Fs13643-016-0263-z/MediaObjects/13643_2016_263_MOESM2_ESM.xlsx",
    }
}

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"

_CLASS_NAMES = ["included", "excluded"]


class SwiftDataset(datasets.GeneratorBasedBuilder):
    """Four datasets (Additional file 1) were generated by the National Toxicology Program (NTP)
    Office of Health Assessment and Translation (OHAT), one dataset (Additional file 2) was provided by
    the Edinburgh CAMARADES group (www.camarades.info)."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    BUILDER_CONFIGS = []
    reviews = REVIEWS
    for dataset_version in reviews:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"swift_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"swift {dataset_version} source schema",
                schema="source",
                subset_id=f"swift_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"swift_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"swift {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"swift_{dataset_version}",
            )
        )

    # Add an "all" config that combines all the reviews -- only for source schema
    BUILDER_CONFIGS.append(
        BigBioConfig(
            name="swift_all_source",
            version=SOURCE_VERSION,
            description=f"swift all source schema",
            schema="source",
            subset_id=f"swift_all",
        )
    )

    DEFAULT_CONFIG_NAME = "swift_all_source"

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

        urls = _URLS[_DATASETNAME]
        data_dir = dl_manager.download(urls)
        pubmed_output_dir = "/".join(self.cache_dir.split("/")[:-3])

        prepare_dataset(input_files=data_dir, output_folder=pubmed_output_dir)

        return [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                gen_kwargs={
                    "split": "train",
                },
            ),
        ]

    def _generate_examples(self, split: str) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""

        data_dir = "/".join(self.cache_dir.split("/")[:-3])
        review = "_".join(self.config.subset_id.split("_")[1:])

        uid = 0

        if review == "all":
            df = pd.DataFrame()
            for r in REVIEWS:
                review_df = pd.read_csv(os.path.join(data_dir, f"{r}.csv"))
                review_df["Review"] = r
                df = pd.concat([df, review_df])
        else:
            df = pd.read_csv(os.path.join(data_dir, f"{review}.csv"))
            df["Review"] = review

        for key, example in df.iterrows():
            review_name = example["Review"]
            title = example["Title"]
            abstract = example["Abstract"]
            label = example["Label"]
            try:
                pmid = str(example["PMID"])
            except (ValueError, KeyError):
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
    x = datasets.load_dataset(__file__, name="swift_BPA_bigbio_text")
    print(x)

    y = datasets.load_dataset(__file__, name="swift_all_source")
    print(y)
