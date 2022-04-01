# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import json
from types import StringTypes
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import defaultdict


client = Config().get_mongo_client(mtype="fungigenome")
db = client[Config().get_mongo_dbname("fungigenome")]


def export_gene_by_geneid(data, option_name, dir_path, bind_obj=None):
    """
    根据gene_id/样品名生成基因序列文件
    :return:
    """
    file_path = os.path.join(dir_path, "%s.fasta" % option_name)
    collection_main = db['gene_predict']
    main = collection_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
    main_id = main['_id']
    collection_seq = db['gene_predict_seq']
    doc = collection_seq.find_one(
        {'predict_id': main_id, 'specimen_id': bind_obj.sheet.option("specimen_id"), 'gene_id': data})
    with open(file_path, 'w') as f:
        f.write(doc['seq_info'] + '\n' + doc['seq_fnn'] + '\n')
    return file_path

