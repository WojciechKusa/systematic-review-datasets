import copy
import xml.etree.ElementTree as ET
from typing import Any, Optional

import pandas as pd


def get_attribute(element: Optional[ET.Element], path: str, default: str) -> str:
    if element is None:
        return default
    result = element.get(path)
    return default if result is None or len(result) == 0 else result


def get_child_texts(
    element: Optional[ET.Element], path: str, default: str, separator: str = " "
) -> str:
    if element is None:
        return default
    result = element.findall(path)

    if result is None or len(result) == 0:
        return default

    all_text: list[str] = []
    for item in result:
        all_text.extend(item.itertext())
    return separator.join(all_text).replace("\n", "")


def extract_review_metadata(root: ET.Element, default: str) -> dict[str, Any]:
    review_metadata: dict[str, Any] = {
        "IS_QUADAS_2": get_attribute(root, "QUADAS2", default=default),
        "REVIEW_TYPE": get_attribute(root, "TYPE", default=default),
        "REVIEW_ID": get_attribute(root, "ID", default=default),
        "REVIEW_DOI": get_attribute(root, "DOI", default=default),
        "REVIEW_GROUP": get_attribute(root, "GROUP_ID", default=default),
    }

    year = get_attribute(
        root.find(".//COVER_SHEET//DATES//LAST_SEARCH//DATE"), "YEAR", default=default
    )
    month = get_attribute(
        root.find(".//COVER_SHEET//DATES//LAST_SEARCH//DATE"), "MONTH", default=default
    )
    review_metadata["LAST_SEARCH"] = f"{month}/{year} last searched"

    review_metadata["REVIEW_PUBLISHED_YEAR"] = get_attribute(
        root.find(".//COVER_SHEET//DATES//LAST_CITATION_ISSUE"), "YEAR", default=default
    )

    review_metadata["REVIEW_TITLE"] = get_child_texts(
        root.find(".//COVER_SHEET"), "TITLE", default=default
    )

    _included_studies = root.findall(
        ".//STUDIES_AND_REFERENCES//STUDIES//INCLUDED_STUDIES//STUDY"
    )
    review_metadata["NUM_INCLUDED_STUDIES"] = len(_included_studies)
    review_metadata["INCLUDED_STUDIES_IDS"] = ";".join(
        [get_attribute(elem, "ID", default=default) for elem in _included_studies]
    )

    return review_metadata


def extract_included_studies(root: ET.Element, default: str) -> list[dict[str, Any]]:
    _included_studies = []

    for elem in root.findall(
        ".//STUDIES_AND_REFERENCES//STUDIES//INCLUDED_STUDIES//STUDY"
    ):
        try:
            _year = int(get_attribute(elem, "YEAR", default=default))
        except ValueError:
            if len(get_attribute(elem, "ID", default=default).split("-")) == 3:
                _year = int(
                    get_attribute(elem, "ID", default=default).split("-")[-1][:4]
                )
            elif len(get_attribute(elem, "ID", default=default).split("-")) == 4:
                _year = int(
                    get_attribute(elem, "ID", default=default).split("-")[-2][:4]
                )
            else:
                _year = 0

        _included_studies.append(
            {
                "STUDY_ID": get_attribute(elem, "ID", default=default),
                "STUDY_PUBLISHED_YEAR": _year,
                "STUDY_NAME": get_attribute(elem, "NAME", default=default),
                "STUDY_DATA_SOURCE": get_attribute(
                    elem, "DATA_SOURCE", default=default
                ),
                "STUDY_DOI": get_attribute(elem, "DOI", default=default),
            }
        )
    return _included_studies


def extract_bias(root: ET.Element, review_data, default: str) -> list[dict[str, Any]]:
    bias_items = root.findall(".//QUALITY_ITEMS//QUALITY_ITEM")
    all_data = []

    for elem in bias_items:
        basic = {
            "BIAS_ELEMENT_ID": get_attribute(elem, "ID", default=default),
            "BIAS_LEVEL_OF_ASSESSMENT": get_attribute(elem, "LEVEL", default=default),
            "CHARACTERISTIC": get_attribute(elem, "CHARACTERISTIC", default=default),
            "IS_CORE_ITEM": get_attribute(elem, "CORE_ITEM", default=default),
            "DOMAIN_NR": get_attribute(elem, "DOMAIN", default=default),
            "DOMAIN_NAME": get_attribute(elem, "DOMAIN_NAME", default=default),
            "SIGNALLING_QUESTION": get_child_texts(elem, "NAME", default=default),
            "SIGNALLING_QUESTION_DESCRIPTION": get_child_texts(
                elem, "DESCRIPTION", separator=" ", default=default
            ),
        }

        bias_data = elem.findall(".//QUALITY_ITEM_DATA")
        for d in bias_data:
            bias_entries = d.findall(".//QUALITY_ITEM_DATA_ENTRY")
            for entry in bias_entries:
                result = copy.deepcopy(basic)
                result["RESULT"] = get_attribute(entry, "RESULT", default=default)
                result["JUDGEMENT_TEXT"] = get_child_texts(
                    entry, "DESCRIPTION", separator=" ", default=default
                )

                result["STUDY_ID"] = get_attribute(entry, "STUDY_ID", default=default)

                result["REVIEW_TITLE"] = review_data["REVIEW_TITLE"]
                result["REVIEW_PUBLISHED_YEAR"] = review_data["REVIEW_PUBLISHED_YEAR"]
                result["REVIEW_DOI"] = review_data["REVIEW_DOI"]
                result["IS_QUADAS_2"] = review_data["IS_QUADAS_2"]

                all_data.append(result)
    return all_data


def get_common_statistics(entry, default, review_data) -> dict[str, str]:
    extracted = {
        "REVIEW_TITLE": review_data["REVIEW_TITLE"],
        "REVIEW_ID": review_data["REVIEW_ID"],
        "REVIEW_PUBLISHED_YEAR": review_data["REVIEW_PUBLISHED_YEAR"],
        "REVIEW_DOI": review_data["REVIEW_DOI"],
        "STUDY_ID": get_attribute(entry, "STUDY_ID", default=default),
    }
    extracted["STUDY_PUBLISHED_YEAR"] = extracted["STUDY_ID"].split("-")[-1]

    return extracted


def extract_statistics_diagnostic(
    root, review_data, default: str
) -> list[dict[str, Any]]:
    stats_items = root.findall(".//ANALYSES_AND_DATA//TESTS//TEST")
    all_data = []

    for elem in stats_items:
        basic = {
            "STATS_ELEMENT_ID": get_attribute(elem, "ID", default=default),
            "SUB_NAME": get_child_texts(elem, "NAME", default=default),
            "FULL_NAME": get_child_texts(elem, "FULL_NAME", default=default),
            "DESCRIPTION": get_child_texts(
                elem, "DESCRIPTION", separator=" ", default=default
            ),
            "TEST_ID": get_child_texts(elem, "ID", default=default),
        }

        test_data = elem.findall(".//TEST_DATA")
        for d in test_data:
            test_entries = d.findall(".//TEST_DATA_ENTRY")
            for entry in test_entries:
                result = copy.deepcopy(basic)

                result["FN"] = get_attribute(entry, "FN", default=default)
                result["FP"] = get_attribute(entry, "FP", default=default)
                result["TN"] = get_attribute(entry, "TN", default=default)
                result["TP"] = get_attribute(entry, "TP", default=default)
                result["SENSITIVITY"] = get_attribute(
                    entry, "SENSITIVITY", default=default
                )
                result["SENS_CI_START"] = get_attribute(
                    entry, "SENS_CI_START", default=default
                )
                result["SENS_CI_END"] = get_attribute(
                    entry, "SENS_CI_END", default=default
                )
                result["SPECIFICITY"] = get_attribute(
                    entry, "SPECIFICITY", default=default
                )
                result["SPEC_CI_START"] = get_attribute(
                    entry, "SPEC_CI_START", default=default
                )
                result["SPEC_CI_END"] = get_attribute(
                    entry, "SPEC_CI_END", default=default
                )

                result |= get_common_statistics(
                    entry=entry, default=default, review_data=review_data
                )

                all_data.append(result)

    return all_data


def extract_study_outcomes(study_root: ET.Element) -> dict[str, str]:
    return dict(study_root.items())


def extract_subgroup_studies(
    subgroup_root: ET.Element, data_tag: str, default: str
) -> list[dict[str, Any]]:
    subgroups_stats = []
    subgroup_stats = {f"SUBGROUP_{key}": value for key, value in subgroup_root.items()}
    subgroup_stats["SUBGROUP_ID"] = subgroup_root.attrib["NO"]
    subgroup_stats["SUBGROUP_NAME"] = get_child_texts(
        subgroup_root, "NAME", default=default
    )
    for study in subgroup_root.findall(data_tag):
        study_outcomes = extract_study_outcomes(study_root=study)
        study_outcomes |= subgroup_stats
        subgroups_stats.append(study_outcomes)
    return subgroups_stats


def extract_outcomes_from_comparison(
    comparison_root: ET.Element, default: str
) -> list[dict[str, Any]]:
    outcomes_stats = []
    outcome_types = {
        "dichotomous": {
            "outcome_tag": "DICH_OUTCOME",
            "subgroup_tag": "DICH_SUBGROUP",
            "data_tag": "DICH_DATA",
        },
        "continuous": {
            "outcome_tag": "CONT_OUTCOME",
            "subgroup_tag": "CONT_SUBGROUP",
            "data_tag": "CONT_DATA",
        },
    }
    for tags in outcome_types.values():
        for outcome in comparison_root.findall(tags["outcome_tag"]):
            outcome_stats = {f"OUTCOME_{key}": value for key, value in outcome.items()}

            outcome_stats["OUTCOME_NO"] = get_attribute(outcome, "NO", default=default)
            outcome_stats["OUTCOME_TYPE"] = outcome.tag
            outcome_stats["OUTCOME_NAME"] = get_child_texts(
                outcome, "NAME", default=default
            )
            outcome_stats["GROUP_LABEL_1"] = get_child_texts(
                outcome, "GROUP_LABEL_1", default=default
            )
            outcome_stats["GROUP_LABEL_2"] = get_child_texts(
                outcome, "GROUP_LABEL_2", default=default
            )
            outcome_stats["GRAPH_LABEL_1"] = get_child_texts(
                outcome, "GRAPH_LABEL_1", default=default
            )
            outcome_stats["GRAPH_LABEL_2"] = get_child_texts(
                outcome, "GRAPH_LABEL_2", default=default
            )
            if outcome.get("SUBGROUPS") == "YES":
                for subgroup in outcome.findall(tags["subgroup_tag"]):
                    subgroups_outcomes = extract_subgroup_studies(
                        subgroup_root=subgroup,
                        data_tag=tags["data_tag"],
                        default=default,
                    )
                    for subgroup_outcome in subgroups_outcomes:
                        subgroup_outcome |= outcome_stats
                        outcomes_stats.append(subgroup_outcome)
            else:
                for study in outcome.findall(tags["data_tag"]):
                    study_outcomes = extract_study_outcomes(study_root=study)

                    study_outcomes |= outcome_stats
                    study_outcomes["SUBGROUP_NAME"] = default
                    study_outcomes["SUBGROUP_ID"] = default
                    study_outcomes["SUBGROUP_EFFECT_SIZE"] = default
                    outcomes_stats.append(study_outcomes)

    return outcomes_stats


def extract_comparisons(root: ET.Element) -> list[dict[str, Any]]:
    comparisons = root.findall("COMPARISON")
    all_data = []

    for comparison in comparisons:
        comparison_id = get_attribute(comparison, "ID", default="")
        comparison_no = get_attribute(comparison, "NO", default="")
        comparison_name = get_child_texts(comparison, "NAME", default="")

        outcomes = extract_outcomes_from_comparison(comparison, default="")

        for outcome in outcomes:
            outcome["COMPARISON_ID"] = comparison_id
            outcome["COMPARISON_NO"] = comparison_no
            outcome["COMPARISON_NAME"] = comparison_name
            all_data.append(outcome)

    return all_data


def extract_statistics_intervention(
    root: ET.Element, review_data: dict[str, str], default: str
) -> list[dict[str, Any]]:
    analyses = root.find(".//ANALYSES_AND_DATA")
    if get_attribute(analyses, "CALCULATED_DATA", default="NO") == "NO":
        return []
    else:
        return extract_comparisons(analyses)


def extract_statistics(
    root: ET.Element, review_data: dict[str, str], default: str
) -> list[dict[str, Any]]:
    if review_data["REVIEW_TYPE"] == "DIAGNOSTIC":
        return extract_statistics_diagnostic(
            root=root, review_data=review_data, default=default
        )

    elif review_data["REVIEW_TYPE"] == "INTERVENTION":
        return extract_statistics_intervention(
            root=root, review_data=review_data, default=default
        )

    else:
        return []


def parse_revman_file(
    xml_path: str,
) -> tuple[
    dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]
]:
    """Parse a RevMan MD5 file and extract the data."""
    _default_value = ""
    root = ET.parse(xml_path).getroot()
    review_data = extract_review_metadata(root, default=_default_value)
    included_studies = extract_included_studies(root, default=_default_value)
    bias_items = extract_bias(root, review_data, default=_default_value)
    stats_items = extract_statistics(root, review_data, default=_default_value)

    return review_data, included_studies, bias_items, stats_items


def save_output(items: list[dict[str, Any]], output_path: str) -> None:
    """Save the output to a CSV file."""
    _column_names = set()
    for _item in items:
        _column_names.update(list(_item.keys()))

    _out_df = pd.DataFrame(items, columns=list(_column_names))
    _out_df.to_csv(output_path, index=False)
