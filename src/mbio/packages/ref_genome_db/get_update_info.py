# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from optparse import OptionParser
import os
from biocluster.config import Config
from bson.objectid import ObjectId
import types
import json
import re
from types import StringTypes
import gridfs
from collections import OrderedDict
import pandas as pd
from biocluster.file import getsize, exists
from biocluster.file import download
import shutil
from biocluster.api.file.lib.transfer import MultiFileTransfer
import sys

project_type = 'ref_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


parser = OptionParser(description='get updata info')
parser.add_option('-i', '--genome_id', dest='genome_id' ,help='start genome_id in sg_genome_db')
parser.add_option('-o', '--output_prefix', dest='output_prefix', help='Output file prefix')
(opts, args) = parser.parse_args()

def get_update_list(genome_id, output_prefix):
    if genome_id.startswith("GM"):
        genome_id = int(genome_id[2:])
    else:
        genome_id = int(genome_id)
    collection = db["sg_genome_db"]
    documents = collection.find()
    output_cmd = output_prefix + "_dir_path.txt"
    output_info = output_prefix + "_info.xls"
    with open(output_cmd, "w") as w1, open(output_info, "w") as w2:
        for document in documents:
            genome_id_1 = int(document["genome_id"][2:])
            if genome_id_1 > genome_id:
                path = "~/app/database/Genome_DB_finish/" + document["anno_path_v2"][:-14]
                w1.write("~/sg-users/liubinxu/script/syn2nb2.sh " + path + "\n")
                common_name = document["common_name"]
                organism_name = document["organism_name"]
                assembly = document["assembly"]
                annot_version = document["annot_version"]
                level = document["level"]
                ensemble_web = document["ensemble_web"]
                try:
                    document["lnc_dir"]
                except:
                    lnc_dir = "否"
                else:
                    lnc_dir = "是"
                if "ncbi" in ensemble_web:
                    source = "NCBI"
                elif "ensembl" in ensemble_web:
                    source = "Ensembl"
                else:
                    source = "other"
                w2.write("\t".join([common_name, organism_name, ensemble_web, assembly, annot_version, source, level, lnc_dir]) + "\n")


if __name__ == '__main__':
    if opts.genome_id and opts.output_prefix:
        get_update_list(opts.genome_id, opts.output_prefix)
    else:
        parser.print_help()