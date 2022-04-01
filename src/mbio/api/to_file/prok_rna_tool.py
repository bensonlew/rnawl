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
import unicodedata
from biocluster.file import getsize, exists
from biocluster.file import download

project_type = 'prok_rna'
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

def export_network_diff(data, option_name, dir_path, bind_obj=None):
    exp_id = bind_obj.sheet.option('exp_id')
    network_path = os.path.join(dir_path, "network.txt")
    diff_seq_path = os.path.join(dir_path, "diff_seq.txt")
    conn = db['sg_diff_detail']
    if not isinstance(data, ObjectId):
        data = ObjectId(data)
    cmp_list = bind_obj.sheet.option('cmp_list')
    type = bind_obj.sheet.option('gene_type')
    significant = bind_obj.sheet.option('significant')
    regulate = bind_obj.sheet.option('regulate')
    pvalue_padjust = bind_obj.sheet.option('stat_type')
    stat_type_operation = bind_obj.sheet.option('stat_type_operation')
    stat_type_value = bind_obj.sheet.option('stat_type_value')
    fc_log2fc = bind_obj.sheet.option('log2fc')
    log2fc_operation = bind_obj.sheet.option('log2fc_operation')
    log2fc_value = bind_obj.sheet.option('log2fc_value')
    geneset_id = bind_obj.sheet.option('geneset_id')
    find_params = dict()
    find_params['diff_id'] = data
    if type == 'mRNA+sRNA':
        pass
    else:
        find_params['type'] = str(type)
    find_params['compare'] = str(cmp_list)
    if significant != "" and not None:
        find_params['significant'] = str(significant)
    if regulate != "" and not None:
        find_params['regulate'] = str(regulate)
    if stat_type_value != "" and not None:
        find_params[pvalue_padjust] = {stat_type_operation: float(stat_type_value)}
    if log2fc_value != "" and not None:
        find_params[fc_log2fc] = {log2fc_operation: float(log2fc_value)}
    if geneset_id != "" and not None:
        conn_ = db['sg_geneset_detail']
        record_ = conn_.find_one({'geneset_id': ObjectId(geneset_id)})
        seq_list = record_['seq_list']
        find_params['seq_id'] = {"$in": seq_list}
    bind_obj.logger.debug("{}".format(find_params))
    print find_params
    record = conn.find(find_params, {'_id': 0, 'seq_id': 1})
    diff_seq = pd.DataFrame(list(record))
    diff_seq.to_csv(diff_seq_path, header=True, index=False, sep='\t')
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    print group_dict
    samples = list()
    cmp_group = cmp_list.strip().split('|')
    for each in group_dict:
        if each in cmp_group:
            samples += group_dict[each]
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    if not isinstance(exp_id, ObjectId):
        exp_id = ObjectId(exp_id)
    conn_exp = db['sg_exp_detail']
    exp_records = conn_exp.find({"exp_id": exp_id},  target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix_filter = exp_matrix[exp_matrix['seq_id'].isin(diff_seq['seq_id'].tolist())]
    exp_matrix_filter.set_index('seq_id', inplace=True)
    exp_matrix_filter = exp_matrix_filter.round(0)
    exp_matrix_filter = exp_matrix_filter.astype(int)
    exp_matrix_filter.reset_index(inplace=True)
    exp_matrix_filter.to_csv(network_path, header=True, index=False, sep='\t')
    return network_path


def export_geneset(data, option_name, dir_path, bind_obj=None):
    geneset_id = ObjectId(data)
    task_id = bind_obj.sheet.option('task_id')
    file_path = os.path.join(dir_path, "geneset_file.txt")
    bind_obj.logger.info('正在导出表达量表')

    conn_ = db['sg_geneset_detail']
    record_ = conn_.find_one({'geneset_id': geneset_id})
    seq_list = record_['seq_list']

    annot_main = db['sg_annotation_query']
    annot_id = annot_main.find_one({"task_id": task_id})['main_id']
    target_cols = OrderedDict(gene_id=1, gene_name=1, _id=0)
    conn_annot = db['sg_annotation_query_detail']
    annot_records = conn_annot.find({"query_id": annot_id}, target_cols)
    annot_matrix = pd.DataFrame(list(annot_records))
    annot_matrix_filter = annot_matrix[annot_matrix['gene_id'].isin(seq_list)]
    annot_matrix_filter.to_csv(file_path, header=True, index=False, sep='\t')
    return file_path


def export_go_circ_exp(data, option_name, dir_path, bind_obj=None):
    diff_id, cmp_group = data.split(';')
    file_path = os.path.join(dir_path, "diff_exp.txt")
    bind_obj.logger.info('正在导出差异结果表')

    target_cols = OrderedDict(seq_id=1, log2fc=1, fc=1, compare=1, _id=0)
    conn_diff = db['sg_diff_detail']
    diff_records = conn_diff.find({'diff_id': ObjectId(diff_id)}, target_cols)
    diff_matrix = pd.DataFrame(list(diff_records))
    diff_matrix_filter = diff_matrix[diff_matrix['compare'] == cmp_group]
    diff_matrix_filter.drop(columns=['compare'], inplace=True)
    diff_matrix_filter['gene_name'] = diff_matrix_filter['seq_id']
    diff_matrix_filter.set_index('seq_id', inplace=True)
    diff_matrix_filter.to_csv(file_path, header=True, index=True, sep='\t')
    return file_path


def export_go_circ_enrich(data, option_name, dir_path, bind_obj=None):
    enrich_id = ObjectId(data)
    file_path = os.path.join(dir_path, "go_file.txt")
    bind_obj.logger.info('正在导出富集结果表')

    target_cols = OrderedDict(go_id=1, go_type=1, description=1, p_corrected=1, seq_list=1, _id=0)
    conn_go = db['sg_geneset_go_enrich_detail']
    go_records = conn_go.find({'go_enrich_id': enrich_id}, target_cols)
    go_matrix = pd.DataFrame(list(go_records))
    go_matrix.rename(columns={'description': 'discription'}, inplace=True)
    go_matrix.set_index('go_id', inplace=True)
    go_matrix.to_csv(file_path, header=True, index=True, sep='\t')
    return file_path


def export_wgcna_exp(data, option_name, dir_path, bind_obj=None):
    exp_id, geneset_id = data.split(';')
    group_dict = bind_obj.sheet.option('group_dict')
    exp_path = os.path.join(dir_path, "wgcna_exp.txt")
    bind_obj.logger.info('正在导出表达量表')

    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]

    conn_ = db['sg_geneset_detail']
    record_ = conn_.find_one({'geneset_id': ObjectId(geneset_id)})
    seq_list = record_['seq_list']

    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    conn_exp = db['sg_exp_detail']
    exp_records = conn_exp.find({"exp_id": ObjectId(exp_id)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix_filter = exp_matrix[exp_matrix['seq_id'].isin(seq_list)]
    exp_matrix_filter.set_index('seq_id', inplace=True)
    exp_matrix_filter.to_csv(exp_path, header=True, index=True, sep='\t')
    return exp_path

