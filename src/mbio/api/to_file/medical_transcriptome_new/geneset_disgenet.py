# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from __future__ import division
import os
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import OrderedDict
import types
import sys
import pandas as pd
import json


project_type = 'medical_transcriptome'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


def export_gene_list_entrez(data, option_name, dir_path, bind_obj=None):
    task_id, geneset_id = data.split(";")
    bind_obj.logger.debug("正在导出基因集")
    gene_list_path = os.path.join(dir_path, "{}_entrez.list".format(geneset_id))
    exp_collection = db["sg_exp"]
    detail_collection = db["sg_exp_detail"]
    geneset_collection = db["sg_geneset_detail"]
    if not isinstance(geneset_id, ObjectId):
        if isinstance(geneset_id, types.StringTypes):
            geneset_id = ObjectId(geneset_id)
        else:
            bind_obj.set_error('geneset_id必须为ObjectId对象或其对应的字符串!', )
    geneset = geneset_collection.find_one({"geneset_id": geneset_id})
    gene_list = geneset["seq_list"]
    exp_main = exp_collection.find_one({"task_id": task_id, "level": "G", "is_rmbe": bool(""), })
    if not exp_main:
        bind_obj.set_error("意外错误，task_id:%s 基因集在sg_exp中未找到！", variables=(geneset_id,), )
    exp_id = exp_main["main_id"]
    entrez_list = []
    for each in gene_list:
        detail = detail_collection.find_one({"exp_id": exp_id, "gene_id": each})
        if not detail:
            continue
        if detail['entrez'] == "":
            continue
        elif detail['entrez'] == "\N":
            continue
        elif detail["entrez"]:
            entrez_list.append(detail["entrez"])
    print(entrez_list)
    with open(gene_list_path, "w") as f:
        for id in entrez_list:
            f.write(id + "\n")
    return gene_list_path


def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} - {}'.format(k, v))
    else:
        pass


def export_entrez_list(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    target_cols = OrderedDict(gene_id=1, _id=0, entrez=1, gene_name=1)
    main = db['sg_exp']
    exp_id = main.find_one({'task_id': data, "level": "G", "is_rmbe": bool("")})['main_id']
    collection = db['sg_exp_detail']
    print target_cols
    cursor = collection.find({'exp_id': exp_id, }, target_cols)
    output = os.path.join(dir_path, 'gene2entrez.txt')
    df = pd.DataFrame(list(cursor))
    exp_matrix = df.set_index('gene_id')
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    # df.to_csv(output, sep='\t', index=False)
    return output



