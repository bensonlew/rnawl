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

id_mapping_file = sys.argv[1]
split_dict = dict()
# os.mkdir("split")

uni_id_last = ""
taxon_id = ""
uni_list = list()

with gzip.open(id_mapping_file, "r") as f:
    for line in f:
        uni_id, db, db_id = line.strip().split("\t")
        if uni_id != uni_id_last and uni_id_last != "":
            if taxon_id in split_dict:
                split_dict[taxon_id].extend(uni_list)
            else:
                split_dict[taxon_id] = uni_list
            uni_list = []

        uni_id_last = uni_id
        uni_list.append(line)
        if db == "NCBI_TaxID":
            taxon_id = db_id
    if taxon_id in split_dict:
        split_dict[taxon_id].extend(uni_list)
    else:
        split_dict[taxon_id] = uni_list

for t in split_dict:
    with gzip.open("split/" +  t + ".id_mapping.gz", "w") as fo:
        fo.write("".join(split_dict[t]))


