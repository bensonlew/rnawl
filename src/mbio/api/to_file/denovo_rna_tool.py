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

project_type = 'denovo_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

def export_upset_venn(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "upset_venn.txt")
    bind_obj.logger.debug("正在导出基因集")
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    geneset_id_list = data.split(',')
    gene_list_all = list()
    gene_name_seq = dict()
    for i in geneset_id_list:
        my_result = main_collection.find_one({'main_id': ObjectId(i)})
        if not my_result:
            bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(i)))
        geneset_name = my_result['name']
        results = collection.find_one({"geneset_id": ObjectId(i)})
        gene_list = results['seq_list']
        gene_name_seq[geneset_name] = gene_list
        gene_list_all.extend(gene_list)
    df = pd.DataFrame.from_dict(gene_name_seq, orient='index').T
    df.to_csv(gene_list_path, header=True, index=False, sep='\t')
    # with open(gene_list_path, 'w') as f:
    #     a = gene_name_seq.keys()
    #     b = '\t'.join(a)
    #     f.write(b + '\n')
    #     for i in gene_list_all:
    #         num_list = list()
    #         for j in a:
    #             if i in gene_name_seq[j]:
    #                 num_list.append(str(1))
    #             else:
    #                 num_list.append("")
    #         f.write('\t'.join(num_list) + '\n')
    return gene_list_path

def export_diff_plot(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析gene_id,log2FC,pvalue
    '''
    diff_id = bind_obj.sheet.option('diff_id')
    compare_group = bind_obj.sheet.option('compare')
    target_cols = OrderedDict(seq_id=1, log2fc=1, pvalue=1, _id=0)
    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id), "compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    columnc_order = ['seq_id', 'log2fc', 'pvalue']
    diff_exp_matrix = diff_exp_matrix[columnc_order].dropna(axis=0)
    output = os.path.join(dir_path, 'diff_plot.txt')
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export diffexpression matrix')
    return output

