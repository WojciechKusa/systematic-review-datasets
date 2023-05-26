import json
import xml.etree.ElementTree as ET
from typing import Tuple

import grobid_tei_xml


def parse_xml_to_dict(element):
    if len(element) == 0:
        return element.text
    result = {}
    for child in element:
        child_data = parse_xml_to_dict(child)
        if child.tag in result:
            if type(result[child.tag]) is list:
                result[child.tag].append(child_data)
            else:
                result[child.tag] = [result[child.tag], child_data]
        else:
            result[child.tag] = child_data
    return result


def convert_xml_to_json(xml_file_path):
    tree = ET.parse(xml_file_path)
    root = tree.getroot()
    data = parse_xml_to_dict(root)
    json_data = json.dumps(data, indent=4)
    return json_data


def get_main_text(xml_tei_path) -> Tuple[str, str, str]:
    with open(xml_tei_path) as f:
        paper = grobid_tei_xml.parse_document_xml(f.read())
        paper = paper.to_dict()

    try:
        title = paper["header"]["title"]
    except KeyError:
        title = ""
    try:
        abstract = paper["abstract"]
    except KeyError:
        abstract = ""
    try:
        main_text = paper["body"]
    except KeyError:
        main_text = ""

    return title, abstract, main_text
