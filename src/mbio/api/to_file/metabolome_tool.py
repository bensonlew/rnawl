# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
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
import unicodedata
from biocluster.file import getsize, exists
from biocluster.file import download

project_type = 'metabolome'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


def export_metabset_exp(data, option_name, dir_path, bind_obj=None):
    exp_id, metabset_id, table_type = data.split(';')
    group_dict = bind_obj.sheet.option('group_dict')
    metab_path = os.path.join(dir_path, "metab_table.txt")
    bind_obj.logger.info('正在导出表达量表')

    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]

    conn_ = db['metab_set_detail']
    record_ = conn_.find_one({'set_id': ObjectId(metabset_id)})
    set_list = record_['set_list']

    target_cols = OrderedDict(metab_id=1, metab=1, _id=0)
    for each in samples:
        target_cols[each] = 1

    exp_detail_table = 'exp_' + table_type
    bind_obj.logger.info(exp_detail_table)
    conn_exp = db[exp_detail_table]
    exp_records = conn_exp.find({"exp_id": ObjectId(exp_id)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix_filter = exp_matrix[exp_matrix['metab_id'].isin(set_list)]
    exp_matrix_filter.set_index('metab', inplace=True)
    exp_matrix_filter.drop(['metab_id'], axis=1, inplace=True)
    exp_matrix_filter.to_csv(metab_path, header=True, index=True, sep='\t')
    return metab_path