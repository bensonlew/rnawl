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

def export_masigpro_geneset(data, option_name, dir_path, bind_obj):
    geneset_id, category, level = data.strip().split(",")
    chk_parm_func(
        sys._getframe().f_code.co_name, data=geneset_id, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    output = os.path.join(dir_path, 'geneset.list')
    if geneset_id in ['All', 'RefAll']:
        collection = db['exp_detail']
        cursor = collection.find({'exp_id': ObjectId(bind_obj.sheet.option('matrix'))})
        if geneset_id == 'All':
            if level == "G":
                dct = {'seq_list': list({document['gene_id'] for document in cursor}), 'category_list': list({document['category'] for document in cursor})}
            else:
                dct = {'seq_list': list({document['transcript_id'] for document in cursor}), 'category_list': list({document['category'] for document in cursor})}
        elif geneset_id == 'RefAll':
            if level == "G":
                dct = {'seq_list': list({document['gene_id'] for document in cursor if document['kind'] == "ref"}), 'category_list': list({document['category'] for document in cursor if document['kind'] == "ref"})}
            else:
                dct = {'seq_list': list({document['transcript_id'] for document in cursor if document['kind'] == "ref"}), 'category_list': list({document['category'] for document in cursor if document['kind'] == "ref"})}
    else:
        collection = db['geneset_detail']
        dct = collection.find_one({'geneset_id': ObjectId(geneset_id)})
    with open(output, "w") as w:
        seq_list = dct['seq_list']
        category_list = dct['category_list']
        for i in range(0, len(category_list)):
            if category_list[i] == category:
                w.writelines('{}\n'.format(seq_list[i]))
    return output

def export_masigpro_matrix(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    collection = db['exp_detail']
    cursor = collection.find({'exp_id': ObjectId(data)})
    output = os.path.join(dir_path, 'matrix.tsv')
    df = pd.DataFrame(list(cursor))
    df.to_csv(output, sep='\t', index=False)
    return output

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
    results = collection.find({"exp_id": ObjectId(query_id), "way": "tpm"})
    output = os.path.join(dir_path, "geneset_info.txt")
    with open(output, "w") as f:
        seq_list = list()
        f.write('seq_id\tcategory\tkind\n')
        for result in results:
            if level == "G":
                seq_id = result["gene_id"]
            else:
                seq_id = result["transcript_id"]
            if seq_id not in seq_list:
                seq_list.append(seq_id)
                category = result["category"]
                kind = result["kind"]
                f.write('{}\t{}\t{}\n'.format(seq_id, category, kind))
    return output


def export_geneset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    info = data.split(";")
    if len(info) == 3:
        exp_id, geneset_id, level = data.split(";")
    if len(info) == 4:
        exp_id, geneset_id, level, rna_type = data.split(";")
    if len(info) == 5:
        exp_id, geneset_id, level, rna_type, way = data.split(";")
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    # export group info
    with open(dir_path + '/group_info.txt', 'w') as f:
        f.write('#sample\tgroup\n')
        group_id = bind_obj.sheet.option('group_id').lower()
        for g in group_dict:
            for s in group_dict[g]:
                if group_id == "all":
                    g = s
                f.write('{}\t{}\n'.format(s, g))
    # decide output columns
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    # get geneset
    if "all" not in geneset_id.lower():
        conn = db['sg_geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        if not geneset_records:
            bind_obj.set_error('geneset not found by query: %s', variables=(geneset_id), code="51008154")
        geneset = geneset_records['seq_list']
    else:
        geneset = list()
    # get all exp matrix
    conn = db['sg_exp_detail']
    # if "refall" in geneset_id.lower():
    #     exp_records = conn.find({"exp_id": ObjectId(exp_id), 'is_new': False}, target_cols)
    # else:
    #     exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    if level == "G":
        exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    if level == "T":
        if rna_type == "all":
            exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
        if rna_type == "mRNA":
            if way == "FPKM":
                exp_records = conn.find({"exp_id": ObjectId(exp_id), "category": "mRNA", "way": "fpkm"}, target_cols)
            if way == "TPM":
                exp_records = conn.find({"exp_id": ObjectId(exp_id), "category": "mRNA", "way": "tpm"}, target_cols)
        if rna_type == "lncRNA":
            if way == "FPKM":
                exp_records = conn.find({"exp_id": ObjectId(exp_id), "category": "lncRNA", "way": "fpkm"}, target_cols)
            if way == "TPM":
                exp_records = conn.find({"exp_id": ObjectId(exp_id), "category": "lncRNA", "way": "tpm"}, target_cols)
        if rna_type == "miRNA":
            exp_records = conn.find({"exp_id": ObjectId(exp_id), "category": "miRNA"}, target_cols)
        if rna_type == "circRNA":
            exp_records = conn.find({"exp_id": ObjectId(exp_id), "category": "circRNA"}, target_cols)

    exp_matrix = pd.DataFrame(list(exp_records))
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: %s', variables=(exp_id), code="51008155")
    exp_matrix = exp_matrix.set_index('seq_id')
    if geneset:
        exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output


def export_group(data, option_name, dir_path, bind_obj=None):
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    group_out = os.path.join(dir_path, option_name)
    with open(group_out, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in group_dict:
            for each in group_dict[key]:
                f.write('{}\t{}\n'.format(each, key))
    return group_out


def get_gene_detail_new(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip().split("|")[0]
    geneset_type = data.strip().split("|")[1]
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", 'status': 'end'})
    if not annot_main:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    if not annot_main:
        bind_obj.set_error("Not Found in sg_annotation_query by query %s", variables=(task_id), code="51008162")
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_annotation_query_detail']
    query_dict = dict(query_id=annot_main_id, )
    if geneset_type.lower() == "t":
        result_dict = dict(_id=0, gene_name=1, gene_id=1, description=1, transcript_id=1)
        result = annot_detail.find(query_dict, result_dict)
        result_pd = pd.DataFrame(list(result))
        result_pd.set_index("transcript_id", inplace=True)
        result_pd = result_pd.loc[:, ["gene_id", "gene_name", "description"]]
        result_pd.columns = ["gene_id", "gene_name", "gene_desc"]
        result_pd["gene_name"][result_pd['gene_name'] == "-"] = None
        result_pd = result_pd.fillna(method="pad", axis=1)
    else:
        result_dict = dict(_id=0, gene_name=1, gene_id=1, description=1)
        result = annot_detail.find(query_dict, result_dict)
        result_pd = pd.DataFrame(list(result))
        result_pd.set_index("gene_id", inplace=True)
        result_pd = result_pd.loc[:, ["gene_name", "description"]]
        result_pd.columns = ["gene_name", "gene_desc"]
        result_pd["gene_name"][result_pd['gene_name'] == "-"] = None
        result_pd = result_pd.fillna(method="pad", axis=1)
    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    os.system(
        r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot
