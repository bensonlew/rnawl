# -*- coding: utf-8 -*-
# __author__ = 'sanger'

from __future__ import division
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

project_type = 'whole_transcriptome'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

def export_geneset_from_exp(data, option_name, dir_path, bind_obj=None):
    collection = db['exp_detail']
    main_collection = db['exp']
    print data
    task_id, level = data.strip().split(",")
    my_result = main_collection.find_one({'task_id': task_id, 'level': level, 'status': 'end'})
    if not my_result:
        bind_obj.set_error("意外错误，task_id: {}在sg_exp中未找到!".format(task_id))
    if 'main_id' in my_result:
        query_id = my_result["main_id"]
    else:
        query_id = my_result['_id']
    results = collection.find({"exp_id": ObjectId(query_id)})
    output = os.path.join(dir_path, "geneset_info.txt")
    with open(output, "w") as f:
        seq_list = set()
        f.write('seq_id\tcategory\tkind\n')
        for result in results:
            if level == "G":
                if "gene_id" in result:
                    seq_id=result["gene_id"]
                else:
                    pass
            else:
                if "transcript_id" in result:
                    seq_id=result["transcript_id"]
                else:
                    pass
            if seq_id not in seq_list:
                seq_list.add(seq_id)
                category=result["category"]
                kind=result["kind"]
                f.write('{}\t{}\t{}\n'.format(seq_id, category, kind))
    return output


