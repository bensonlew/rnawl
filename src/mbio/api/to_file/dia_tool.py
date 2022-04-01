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

project_type = 'dia'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

def export_upset_venn(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "upset_venn.txt")
    bind_obj.logger.debug("正在导出基因集")
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    proteinset_id_list = data.split(',')
    gene_list_all = list()
    gene_name_seq = dict()
    for i in proteinset_id_list:
        my_result = main_collection.find_one({'main_id': ObjectId(i)})
        if not my_result:
            bind_obj.set_error("意外错误，proteinset_id:{}在sg_proteinset中未找到！".format(ObjectId(i)))
        geneset_name = my_result['name']
        results = collection.find_one({"proteinset_id": ObjectId(i)})
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
    target_cols = OrderedDict(accession_id=1, log2fc=1, pvalue=1, _id=0)
    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id), "compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    columnc_order = ['accession_id', 'log2fc', 'pvalue']
    diff_exp_matrix = diff_exp_matrix[columnc_order].dropna(axis=0)
    output = os.path.join(dir_path, 'diff_plot.txt')
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export diffexpression matrix')
    return output

def export_diff_ma_plot(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析gene_id,log2FC,pvalue
    '''
    diff_id = bind_obj.sheet.option('diff_id')
    compare_group = bind_obj.sheet.option('compare')
    ctrl ,test = compare_group.split("|")[-1],compare_group.split("|")[0]
    target_cols = OrderedDict(accession_id=1, log2fc=1, pvalue=1,_id=0)
    target_cols[ctrl] = 1
    target_cols[test] = 1
    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id), "compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    columnc_order = ['accession_id',test, ctrl,'log2fc', 'pvalue']
    diff_exp_matrix = diff_exp_matrix[columnc_order]
    output = os.path.join(dir_path, 'diff_ma.txt')
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export diffexpression matrix')
    return output


def export_circularbar(data, option_name, dir_path, bind_obj=None):
    bind_obj.logger.debug("正在导出GO分类信息")
    collection = db['sg_annotation_go_detail']
    main_collection = db['sg_annotation_go']
    go_id = data.split('____')[0]
    level = int(data.split('____')[1])
    my_result = main_collection.find_one({'main_id': ObjectId(go_id)})
    if not my_result:
        bind_obj.set_error("意外错误，main_id:{}在sg_annotation_go中未找到！".format(ObjectId(go_id)))
    target_cols = OrderedDict(goterm_2=1, goterm=1, seq_num=1, _id=0)
    results = collection.find({"go_id": ObjectId(go_id), "level": level}, target_cols)
    go_df = pd.DataFrame(list(results))
    columnc_order = ['goterm_2', 'goterm', 'seq_num']
    go_df = go_df[columnc_order]
    # sort the matrix by seq_num and goterm_2
    go_df[["seq_num"]] = go_df[["seq_num"]].astype(int)
    go_df = go_df.sort_values(by=['goterm','seq_num'], ascending=[True,False])
    go_df = pd.concat([i[1][:10] for i in go_df.groupby(by=["goterm"])])
    go_df = go_df.rename(columns={'goterm_2': 'Go Term_2', 'goterm': 'Go Term','seq_num': 'Protein Num'})
    #
    output = os.path.join(dir_path, "circularbar.txt")
    go_df.to_csv(output, sep='\t', header=True, index=False)
    print('success to export go dataframe matrix')
    return output

def export_diff_radar_plot(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析gene_id,log2FC,pvalue
    '''
    diff_id,proteinset_id = bind_obj.sheet.option('diff_id').split(",")
    conn_proset = db['sg_proteinset']
    proteinset_result = conn_proset.find_one({"main_id": ObjectId(proteinset_id)})
    if not proteinset_result:
        bind_obj.set_error("意外错误:未找到基因集_id为{}的基因集信息".format(proteinset_id))
    collection = db['sg_proteinset_detail']
    results = collection.find_one({"proteinset_id": ObjectId(proteinset_id)})
    proteinset_names = set(results["seq_list"])
    all_proteins = list(proteinset_names)
    compare_group = bind_obj.sheet.option('compare')
    ctrl ,test = compare_group.split("|")[-1],compare_group.split("|")[0]
    target_cols = OrderedDict(accession_id=1, log2fc=1,_id=0)
    target_cols[ctrl] = 1
    target_cols[test] = 1
    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id), "compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    columnc_order = ['accession_id',test, ctrl,'log2fc']
    diff_exp_matrix = diff_exp_matrix[columnc_order]
    diff_exp_matrix = diff_exp_matrix.set_index("accession_id").loc[all_proteins]
    output = os.path.join(dir_path, 'diff_radar.txt')
    diff_exp_matrix.to_csv(output, sep='\t', header=True)
    # diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export diffexpression matrix')
    return output

def export_kegg_enrich_table(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析gene_id,log2FC,pvalue
    '''
    enrich_id = bind_obj.sheet.option('enrich_id')
    pvalue_padjust = bind_obj.sheet.option('pvalue_padjust')
    conn = db['sg_proteinset_kegg_enrich_detail']
    target_cols = OrderedDict(id=1, study_count=1, _id=0)
    if pvalue_padjust == "pvalue":
        target_cols["pvalue"] =1
        columnc_order = ["id","study_count","pvalue"]
    else:
        target_cols["corrected_pvalue"] = 1
        columnc_order = ["id", "study_count", "corrected_pvalue"]
    kegg_enrich_records = conn.find({"kegg_enrich_id": ObjectId(enrich_id)}, target_cols)
    kegg_enrich_matrix = pd.DataFrame(list(kegg_enrich_records))
    kegg_enrich_matrix = kegg_enrich_matrix[columnc_order]
    kegg_enrich_matrix["id"] = kegg_enrich_matrix["id"].apply(lambda  x :"map"+x[-5:])
    output = os.path.join(dir_path, 'kegg_enrich_table.txt')
    kegg_enrich_matrix.to_csv(output, sep='\t', header=True,index=False)
    return output

def export_go_enrich_table(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析gene_id,log2FC,pvalue
    '''
    enrich_id = bind_obj.sheet.option('enrich_id')
    pvalue_padjust = bind_obj.sheet.option('pvalue_padjust')
    conn = db['sg_proteinset_go_enrich_detail']
    target_cols = OrderedDict(go_id=1, study_count=1, _id=0)
    if pvalue_padjust == "pvalue":
        target_cols["p_uncorrected"] =1
        columnc_order = ["go_id","study_count","p_uncorrected"]
    else:
        target_cols["p_corrected"] = 1
        columnc_order = ["go_id", "study_count", "p_corrected"]
    go_enrich_records = conn.find({"go_enrich_id": ObjectId(enrich_id)}, target_cols)
    go_enrich_matrix = pd.DataFrame(list(go_enrich_records))
    go_enrich_matrix = go_enrich_matrix[columnc_order]
    # kegg_enrich_matrix["id"] = kegg_enrich_matrix["id"].apply(lambda  x :"map"+x[-5:])
    output = os.path.join(dir_path, 'go_enrich_table.txt')
    go_enrich_matrix.to_csv(output, sep='\t', header=True,index=False)
    return output


def export_proteinset_exp(data, option_name, dir_path, bind_obj=None):
    exp_id, proteinset_id = data.split(';')
    group_dict = bind_obj.sheet.option('group_dict')
    metab_path = os.path.join(dir_path, "metab_table.txt")
    bind_obj.logger.info('正在导出表达量表')

    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]

    conn_ = db['sg_proteinset_detail']
    record_ = conn_.find_one({'proteinset_id': ObjectId(proteinset_id)})
    seq_list = record_['seq_list']

    target_cols = OrderedDict(accession_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    conn_exp = db['sg_express_detail']
    exp_records = conn_exp.find({"express_id": ObjectId(exp_id)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix_filter = exp_matrix[exp_matrix['accession_id'].isin(seq_list)]
    exp_matrix_filter.set_index('accession_id', inplace=True)
    exp_matrix_filter.to_csv(metab_path, header=True, index=True, sep='\t')
    return metab_path