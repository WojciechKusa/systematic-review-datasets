"""This module contains the class for a systematic review.
It is used to store the data of a systematic review and to convert it to a format that can be used by the machine
learning models."""
from typing import Optional


class SystematicReview:
    """This class represents a systematic review."""

    review_id: str
    title: str
    search_results: list[dict[str, str]]  # list of papers found by the search
    first_screening: list[dict[str, str]]  # usually title and abstract screening

    abstract: Optional[str]
    authors: Optional[list[str]]

    search_queries: Optional[list[str]]

    inclusion_criteria: Optional[list[str]]
    exclusion_criteria: Optional[list[str]]

    seed_studies: Optional[list[dict[str, str]]]  # list of seed studies
    second_screening: Optional[list[dict[str, str]]]  # usually full text screening

    def __init__(
        self,
        review_id: str,
        title: str,
        search_results: list[dict[str, str]],
        first_screening: list[dict[str, str]],
        abstract: Optional[str] = None,
        authors: Optional[list[str]] = None,
        search_queries: Optional[list[str]] = None,
        inclusion_criteria: Optional[list[str]] = None,
        exclusion_criteria: Optional[list[str]] = None,
        seed_studies: Optional[list[dict[str, str]]] = None,
        second_screening: Optional[list[dict[str, str]]] = None,
    ):
        self.review_id = review_id
        self.title = title
        self.search_results = search_results
        self.first_screening = first_screening

        self.abstract = abstract
        self.authors = authors
        self.search_queries = search_queries
        self.inclusion_criteria = inclusion_criteria
        self.exclusion_criteria = exclusion_criteria
        self.seed_studies = seed_studies
        self.second_screening = second_screening
