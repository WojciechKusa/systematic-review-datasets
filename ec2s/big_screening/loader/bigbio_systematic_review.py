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

"""
This template serves as a starting point for contributing a dataset to the BigScience Biomedical repo.

When modifying it for your dataset, look for TODO items that offer specific instructions.

Full documentation on writing dataset loading scripts can be found here:
https://huggingface.co/docs/datasets/add_dataset.html

To create a dataset loading script you will create a class and implement 3 methods:
  * `_info`: Establishes the schema for the dataset, and returns a datasets.DatasetInfo object.
  * `_split_generators`: Downloads and extracts data for each split (e.g. train/val/test) or associate local data with each split.
  * `_generate_examples`: Creates examples from data on disk that conform to each schema defined in `_info`.

TODO: Before submitting your script, delete this doc string and replace it with a description of your dataset.

[bigbio_schema_name] = (kb, pairs, qa, text, t2t, entailment)
"""

import os
from typing import List, Tuple, Dict
import importlib

import datasets
from ec2s.loader.bigbiohub import text_features
from ec2s.loader.bigbiohub import BigBioConfig
from ec2s.loader.bigbiohub import Tasks
from ec2s.loader.systematic_review import SystematicReview


_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

# TODO: Add BibTeX citation
_CITATION = """\
@article{,
  author    = {},
  title     = {},
  journal   = {},
  volume    = {},
  year      = {},
  url       = {},
  doi       = {},
  biburl    = {},
  bibsource = {}
}
"""

_DATASETNAME = "ec2s"
_DISPLAYNAME = "EC2S"

_DESCRIPTION = """\
This dataset covers the general task of citation screening.
"""

# TODO: Add a link to an official homepage for the dataset here (if possible)
_HOMEPAGE = ""

# TODO: Add the licence for the dataset here (if possible)
# Note that this doesn't have to be a common open source license.
# Some datasets have custom licenses. In this case, simply put the full license terms
# into `_LICENSE`
_LICENSE = ""

# TODO: Add links to the urls needed to download your dataset files.
#  For local datasets, this variable can be an empty dictionary.
#  For publicly available datasets you will most likely end up passing these URLs to dl_manager in _split_generators.
#  In most cases the URLs will be the same for the source and bigbio config.
#  However, if you need to access different files for each config you can have multiple entries in this dict.
#  This can be an arbitrarily nested dict/list of URLs (see below in `_split_generators` method)
_URLS = {
    _DATASETNAME: {
        "cohen": "https://dmice.ohsu.edu/cohenaa/epc-ir-data/epc-ir.clean.tsv",
    }
}  #  FIXME: leave most probably empty

CONFIG_PATH = os.path.join(os.path.dirname(__file__), "config.json")

_SUPPORTED_TASKS = [Tasks.TEXT_CLASSIFICATION, Tasks.QUESTION_ANSWERING]

_SOURCE_VERSION = "1.0.0"
_BIGBIO_VERSION = "1.0.0"


_CLASS_NAMES = ["included", "excluded"]


# TODO: Name the dataset class to match the script name using CamelCase instead of snake_case
#  Append "Dataset" to the class name: BioASQ --> BioasqDataset
class NewDataset(datasets.GeneratorBasedBuilder):
    """TODO: Short description of my dataset."""

    SOURCE_VERSION = datasets.Version(_SOURCE_VERSION)
    BIGBIO_VERSION = datasets.Version(_BIGBIO_VERSION)

    # You will be able to load the "source" or "bigbio" configurations with
    # ds_source = datasets.load_dataset('my_dataset', name='source')
    # ds_bigbio = datasets.load_dataset('my_dataset', name='bigbio')

    # For local datasets you can make use of the `data_dir` and `data_files` kwargs
    # https://huggingface.co/docs/datasets/add_dataset.html#downloading-data-files-and-organizing-splits
    # ds_source = datasets.load_dataset('my_dataset', name='source', data_dir="/path/to/data/files")
    # ds_bigbio = datasets.load_dataset('my_dataset', name='bigbio', data_dir="/path/to/data/files")

    # TODO: For each dataset, implement Config for Source and BigBio;
    #  If dataset contains more than one subset (see examples/bioasq.py) implement for EACH of them.
    #  Each of them should contain:
    #   - name: should be unique for each dataset config eg. bioasq10b_(source|bigbio)_[bigbio_schema_name]
    #   - version: option = (SOURCE_VERSION|BIGBIO_VERSION)
    #   - description: one line description for the dataset
    #   - schema: options = (source|bigbio_[bigbio_schema_name])
    #   - subset_id: subset id is the canonical name for the dataset (eg. bioasq10b)
    #  where [bigbio_schema_name] = (kb, pairs, qa, text, t2t, entailment)

    BUILDER_CONFIGS = []
    dataset_versions = ["cohen"]
    for dataset_version in dataset_versions:
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"ec2s_{dataset_version}_source",
                version=SOURCE_VERSION,
                description=f"ec2s {dataset_version} source schema",
                schema="source",
                subset_id=f"ec2s_{dataset_version}",
            )
        )
        BUILDER_CONFIGS.append(
            BigBioConfig(
                name=f"ec2s_{dataset_version}_bigbio_text",
                version=BIGBIO_VERSION,
                description=f"ec2s {dataset_version} BigBio schema",
                schema="bigbio_text",
                subset_id=f"ec2s_{dataset_version}",
            )
        )

    DEFAULT_CONFIG_NAME = "ec2s_source"

    def _info(self) -> datasets.DatasetInfo:

        # Create the source schema; this schema will keep all keys/information/labels as close to the original dataset as possible.

        # You can arbitrarily nest lists and dictionaries.
        # For iterables, use lists over tuples or `datasets.Sequence`

        if self.config.schema == "source":
            features = datasets.Features(
                {
                    "document_id": datasets.Value("string"),
                    "text": datasets.Value("string"),
                    "label": [datasets.ClassLabel(names=_CLASS_NAMES)],
                }
            )
            # # TODO: Create your source schema here
            # raise NotImplementedError()
            #
            # # EX: Arbitrary NER type dataset
            # # features = datasets.Features(
            # #    {
            # #        "doc_id": datasets.Value("string"),
            # #        "text": datasets.Value("string"),
            # #        "entities": [
            # #            {
            # #                "offsets": [datasets.Value("int64")],
            # #                "text": datasets.Value("string"),
            # #                "type": datasets.Value("string"),
            # #                "entity_id": datasets.Value("string"),
            # #            }
            # #        ],
            # #    }
            # # )

        # Choose the appropriate bigbio schema for your task and copy it here. You can find information on the schemas in the CONTRIBUTING guide.

        # In rare cases you may get a dataset that supports multiple tasks requiring multiple schemas. In that case you can define multiple bigbio configs with a bigbio_text format.

        elif self.config.schema == "bigbio_text":
            features = text_features

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
        urls = _URLS[_DATASETNAME][dataset_name]
        data_dir = dl_manager.download_and_extract(urls)

        config_file = f'datasets/{dataset_name}/setup.json'
        _preprocessing_module = f'ec2s.datasets.{dataset_name}.{dataset_name}'
        prepare_dataset = getattr(importlib.import_module(_preprocessing_module), 'prepare_dataset')

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
            datasets.SplitGenerator(
                name=datasets.Split.TEST,
                gen_kwargs={
                    "filepath": os.path.join(data_dir, "test.jsonl"),
                    "split": "test",
                },
            ),
            datasets.SplitGenerator(
                name=datasets.Split.VALIDATION,
                gen_kwargs={
                    "filepath": os.path.join(data_dir, "dev.jsonl"),
                    "split": "dev",
                },
            ),
        ]

    # method parameters are unpacked from `gen_kwargs` as given in `_split_generators`

    # TODO: change the args of this function to match the keys in `gen_kwargs`. You may add any necessary kwargs.

    def _generate_examples(self, filepath, split: str) -> Tuple[int, Dict]:
        """Yields examples as (key, example) tuples."""
        # TODO: This method handles input defined in _split_generators to yield (key, example) tuples from the dataset.

        # The `key` is for legacy reasons (tfds) and is not important in itself, but must be unique for each example.

        # NOTE: For local datasets you will have access to self.config.data_dir and self.config.data_files

        if self.config.schema == "source":
            # TODO: yield (key, example) tuples in the original dataset schema
            for key, example in thing:
                yield key, example

        elif self.config.schema == "bigbio_text":
            # TODO: yield (key, example) tuples in the bigbio schema
            for key, example in thing:
                yield key, example


# This template is based on the following template from the datasets package:
# https://github.com/huggingface/datasets/blob/master/templates/new_dataset_script.py


# This allows you to run your dataloader with `python ec2s.py` during development
# TODO: Remove this before making your PR
if __name__ == "__main__":
    datasets.load_dataset(__file__, name="ec2s_cohen_source")
