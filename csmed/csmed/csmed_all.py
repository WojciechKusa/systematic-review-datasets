import json
from typing import Union

import datasets

from csmed.loader.bigbiohub import Tasks
from csmed.loader.bigbiohub import text_features
from csmed.csmed.csmed_basic import CSMED_BASIC_REVIEWS
from csmed.csmed.csmed_cochrane import CSMED_COCHRANE_REVIEWS


ALL_REVIEWS : dict = CSMED_BASIC_REVIEWS
ALL_REVIEWS.update(CSMED_COCHRANE_REVIEWS)

_LANGUAGES = ["English"]
_PUBMED = True
_LOCAL = False

_CITATION = """\
"""

_DATASETNAME = "csmed_all"
_DISPLAYNAME = "csmed_all"

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


class CSMeDCochrane:
    """296 Cochrane systematic reviews from the biomedical domain.
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

        :param base_path: path to the directory where the CSMeDCochrane file is stored
        """
        path = f"{base_path}/../../../data/external/"
        csmed_all_datasets = {"TRAIN": {}, "EVAL": {}}
        for dataset_name in ALL_REVIEWS:
            with open(f"{path}/{dataset_name}/data_index.json") as f:
                collection_details = json.load(f)

            for review_name in ALL_REVIEWS[dataset_name]:
                details_file = (
                    f"{path}/{dataset_name}/{review_name}/review_details.json"
                )
                with open(details_file) as f:
                    dataset_details = json.load(f)

                _dataset = datasets.load_dataset(
                    path=f"{base_path}/../{dataset_name}/{dataset_name}.py",
                    name=f"{dataset_name}_{review_name}_bigbio_text",
                )
                csmed_all_datasets["TRAIN"][review_name] = {
                    "review_name": review_name,
                    "dataset_details": dataset_details,
                    "data": _dataset,
                }

        return csmed_all_datasets

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
