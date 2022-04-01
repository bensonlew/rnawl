# -*- coding: utf-8 -*-
import os
import math
import subprocess
import numpy as np
import pandas as pd
import argparse
from multiprocessing import Pool
import glob
import matplotlib
import gzip
import json
import copy
from mbio.api.database.gene_db import genome
import sys
import xml.etree.ElementTree as ET

rdf_file = sys.argv[1]

a = open(rdf_file, "r")
xml = ET.parse(a)
root = xml.getroot()

rdf_ns = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
des =  root.iterfind(rdf_ns + 'Description')
namespace_dict = {
    "xmlns" : "http://purl.uniprot.org/core/",
    "xmlns:dcterms" : "http://purl.org/dc/terms/",
    "xmlns:rdf" : "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "xmlns:rdfs" : "http://www.w3.org/2000/01/rdf-schema#",
    "xmlns:owl" : "http://www.w3.org/2002/07/owl#",
    "xmlns:skos" : "http://www.w3.org/2004/02/skos/core#",
    "xmlns:bibo" : "http://purl.org/ontology/bibo/",
    "xmlns:foaf" : "http://xmlns.com/foaf/0.1/",
    "xmlns:void" : "http://rdfs.org/ns/void#",
    "xmlns:sd" : "http://www.w3.org/ns/sparql-service-description#",
    "xmlns:faldo" : "http://biohackathon.org/resource/faldo#"
}

all_e = ["Description", "type", "abbreviation", "identifier", "category", "label", "linkIsExplicit", "seeAlso", "urlTemplate"]
print "\t".join(all_e)

for ades in des:
    about = ades.attrib.get(rdf_ns + 'about', '')
    eles = dict()
    
    for chi in ades.iter():
        tag = chi.tag.split("}")[1]
        text = chi.text
        eles[tag] = str(text).strip()
        
    # print eles

    print "\t".join([eles.get(e, "") for e in all_e])

