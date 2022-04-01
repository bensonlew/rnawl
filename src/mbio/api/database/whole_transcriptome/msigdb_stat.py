# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import sys
import os
import shutil
import pandas as pd
from collections import Counter
from collections import OrderedDict
from api_base import ApiBase
from bson.son import SON
import csv

## /mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/all.symbols.gmt
msig_file = sys.argv[1]
msig_df = pd.read_table(msig_file, sep="\t", header=0)
msig_df = msig_df.fillna("")

json_list = list()

for species in ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus']:
    msig_df_spe = msig_df[msig_df["organism"]==species]
    species_json = {"organism": species}
    c1_num = OrderedDict(msig_df_spe.groupby(['c1']).count()['name'])
    class_dict = OrderedDict()
    for c1 in c1_num:
        c2_c3_dict = dict()
        msig_df_spe_c1 = msig_df_spe[msig_df_spe["c1"]==c1]
        c2_c3 = dict(msig_df_spe_c1.groupby(['c2', 'c3']).count()['name'])
        for key,value in c2_c3.items():
            c2_c3_dict["|".join(key)] = {"name": "|".join(key), "num": value}
        class_dict.update({
            c1: {
                "name": c1,
                "num": c1_num[c1],
                "sub_class": c2_c3_dict
            }
        })
    json_list.append(SON({"organism": species, "class": class_dict}))

my_api = ApiBase(None)
my_collection = my_api.db["msigdb_stat"]
my_collection.insert_many(json_list)

with open(msig_file) as in_handler:
    lst = [dic for dic in csv.DictReader(in_handler, delimiter='\t')]
my_api.create_db_table('msigdb', lst)

print json_list
