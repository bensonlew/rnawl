# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET
import xmltodict
import json
import gzip

def parse_all_xml(xml_path, taxon_list = []):
    taxon_dict = dict()
    for taxon in taxon_list:
        taxon_dict[taxon] = dict()

    sub_lines = ""
    in_sub = False
    with gzip.open(xml_path, 'rb') as f:
        for line in f:
            if line.startswith("<entry "):
                sub_lines = line
                in_sub = True
            elif in_sub:
                sub_lines += line

            if line.startswith("</entry>"):
                in_sub = False
                get_sub_xml_info(sub_lines, taxon_list, taxon_dict)

    for k, v in taxon_dict.items():
        with open(k + ".unprot_sprot.json", 'w') as f:
            f.write(json.dumps(v, indent=4, ensure_ascii=False).encode('utf8'))

def get_sub_xml_info(xml_string, taxon_list = [], taxon_dict = dict()):

    data_dict = xmltodict.parse(xml_string)
    entry =  data_dict["entry"]
    try:
        if 'organism' in entry and \
           'dbReference' in entry['organism'] and \
           '@id' in entry['organism']['dbReference']:
            if entry['organism']['dbReference']['@id'] in taxon_list:
                if type(entry["accession"]) == unicode:
                    taxon_dict[entry['organism']['dbReference']['@id']].update(
                        {entry["accession"]: entry}
                    )
                elif type(entry["accession"]) == list:
                    taxon_dict[entry['organism']['dbReference']['@id']].update(
                        {"|".join(entry["accession"]): entry}
                    )
        else:
            print "no taxon id {}".format(entry["accession"])
    except:
        print entry


if __name__ == "__main__":
    xml_path = sys.argv[1]
    taxon_list = sys.argv[2].split(",")
    parse_all_xml(xml_path, taxon_list)
