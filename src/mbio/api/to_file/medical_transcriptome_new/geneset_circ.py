# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from __future__ import division
import os
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import OrderedDict
import pandas as pd
import types
import json


project_type = 'medical_transcriptome'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

def export_go_enrich_matrix(data, option_name, dir_path, bind_obj=None):

    go_enrich_id = bind_obj.sheet.option('enrich_id')
    go_type = bind_obj.sheet.option('go_type')

    target_cols = OrderedDict(go_id=1, go_type=1, discription=1, p_corrected=1, p_uncorrected=1, seq_str=1, depth=1, _id=0)
    bind_obj.logger.debug("导出GO富集表")
    # get geneset
    conn = db['sg_geneset_go_enrich_detail']
    if go_type == 'ALL' or go_type == 'All':
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id)}, target_cols)
    else:
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id), "go_type":go_type}, target_cols)
    go_enrich_matrix = pd.DataFrame(list(go_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in go_enrich_matrix['seq_str']]
    go_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    go_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export go_enrich matrix')
    return output

def export_kegg_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出kegg富集表格
    '''
    kegg_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, id=1, term=1, pvalue=1, corrected_pvalue=1, seq_str=1, kegg_type=1)
    bind_obj.logger.debug("导出KEGG富集表")
    # get geneset
    conn = db['sg_geneset_kegg_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(kegg_enrich_id))
    kegg_enrich_records = conn.find({"kegg_enrich_id": ObjectId(kegg_enrich_id)}, target_cols)
    kegg_enrich_matrix = pd.DataFrame(list(kegg_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in kegg_enrich_matrix['seq_str']]
    kegg_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    kegg_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export kegg_enrich matrix')
    return output

def export_disgenet_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出disgenet富集表格
    '''
    disgenet_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, disease_id=1, disease_name=1, pvalue=1, padjust=1, seq_str=1)
    bind_obj.logger.debug("导出DisGeNET富集表")
    # get geneset
    conn = db['sg_geneset_disgenet_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(disgenet_enrich_id))
    disgenet_enrich_records = conn.find({"disgenet_enrich_id": ObjectId(disgenet_enrich_id)}, target_cols)
    disgenet_enrich_matrix = pd.DataFrame(list(disgenet_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in disgenet_enrich_matrix['seq_str']]
    disgenet_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    disgenet_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export disgenet_enrich matrix')
    return output

def export_do_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出do富集表格
    '''
    do_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, do_id=1, do_name=1, term_type=1, pvalue=1, padjust=1, seq_str=1)
    bind_obj.logger.debug("导出DO富集表")
    # get geneset
    conn = db['sg_geneset_do_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(do_enrich_id))
    do_enrich_records = conn.find({"do_enrich_id": ObjectId(do_enrich_id)}, target_cols)
    do_enrich_matrix = pd.DataFrame(list(do_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in do_enrich_matrix['seq_str']]
    do_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    do_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export do_enrich matrix')
    return output


def export_reactome_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出reactome富集表格
    '''
    reactome_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, pathway_id=1, description=1, pvalue=1, padjust=1, seq_str=1, category=1)
    bind_obj.logger.debug("导出Reactome富集表")
    # get geneset
    conn = db['sg_geneset_reactome_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(reactome_enrich_id))
    reactome_enrich_records = conn.find({"reactome_enrich_id": ObjectId(reactome_enrich_id)}, target_cols)
    reactome_enrich_matrix = pd.DataFrame(list(reactome_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in reactome_enrich_matrix['seq_str']]
    reactome_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    reactome_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export reactome_enrich matrix')
    return output


def get_gene_detail_new(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip().split("|")[0]
    level = data.strip().split("|")[1]
    annot_table = db['sg_exp']
    annot_main = annot_table.find_one({"task_id": task_id, "level": level, "is_rmbe": bool("")})
    if not annot_main:
        bind_obj.set_error("Not Found in sg_exp by query %s and level %s", variables=(task_id,level))
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_exp_detail']
    query_dict = dict(exp_id=annot_main_id, )
    if level.lower() == "t":
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
    os.system(r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot

def export_compare_exp_fc(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析accession_id 和 fc
    '''

    diff_id = bind_obj.sheet.option('diff_id')
    compare_group = bind_obj.sheet.option('compare_group')
    target_cols = OrderedDict(seq_id=1, log2fc=1, _id=0)

    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id),"compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    output = os.path.join(dir_path, option_name)
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export expression matrix')
    return output

def export_diff_go_enrich_matrix(data, option_name, dir_path, bind_obj=None):

    go_enrich_id = bind_obj.sheet.option('enrich_id')
    go_type = bind_obj.sheet.option('go_type')

    target_cols = OrderedDict(go_id=1, go_type=1, discription=1, p_corrected=1, p_uncorrected=1, seq_str=1, depth=1, _id=0)
    bind_obj.logger.debug("导出GO富集表")
    # get geneset
    conn = db['sg_diff_geneset_go_enrich_detail']
    if go_type == 'ALL' or go_type == 'All':
        go_enrich_records = conn.find({"diff_go_enrich_id": ObjectId(go_enrich_id)}, target_cols)
    else:
        go_enrich_records = conn.find({"diff_go_enrich_id": ObjectId(go_enrich_id), "go_type":go_type}, target_cols)
    go_enrich_matrix = pd.DataFrame(list(go_enrich_records))
    try:
        seq_list2 = [i.replace(';', '|') for i in go_enrich_matrix['seq_str']]
    except:
        seq_list2 = ["|".join(i) for i in go_enrich_matrix['seq_str']]
    go_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    go_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export go_enrich matrix')
    return output

def export_diff_kegg_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出kegg富集表格
    '''
    kegg_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, id=1, term=1, pvalue=1, corrected_pvalue=1, seq_str=1, kegg_type=1)
    bind_obj.logger.debug("导出KEGG富集表")
    # get geneset
    conn = db['sg_diff_geneset_kegg_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(kegg_enrich_id))
    kegg_enrich_records = conn.find({"diff_kegg_enrich_id": ObjectId(kegg_enrich_id)}, target_cols)
    kegg_enrich_matrix = pd.DataFrame(list(kegg_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in kegg_enrich_matrix['seq_str']]
    kegg_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    kegg_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export kegg_enrich matrix')
    return output

def export_diff_disgenet_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出disgenet富集表格
    '''
    disgenet_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, disease_id=1, disease_name=1, pvalue=1, padjust=1, seq_str=1)
    bind_obj.logger.debug("导出DisGeNET富集表")
    # get geneset
    conn = db['sg_diff_geneset_disgenet_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(disgenet_enrich_id))
    disgenet_enrich_records = conn.find({"disgenet_enrich_id": ObjectId(disgenet_enrich_id)}, target_cols)
    disgenet_enrich_matrix = pd.DataFrame(list(disgenet_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in disgenet_enrich_matrix['seq_str']]
    disgenet_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    disgenet_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export disgenet_enrich matrix')
    return output

def export_diff_do_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出do富集表格
    '''
    do_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, do_id=1, do_name=1, term_type=1, pvalue=1, padjust=1, seq_str=1)
    bind_obj.logger.debug("导出DO富集表")
    # get geneset
    conn = db['sg_diff_geneset_do_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(do_enrich_id))
    do_enrich_records = conn.find({"diff_do_enrich_id": ObjectId(do_enrich_id)}, target_cols)
    do_enrich_matrix = pd.DataFrame(list(do_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in do_enrich_matrix['seq_str']]
    do_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    do_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export do_enrich matrix')
    return output


def export_diff_reactome_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出reactome富集表格
    '''
    reactome_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, pathway_id=1, description=1, pvalue=1, padjust=1, seq_str=1, category=1)
    bind_obj.logger.debug("导出Reactome富集表")
    # get geneset
    conn = db['sg_diff_geneset_reactome_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(reactome_enrich_id))
    reactome_enrich_records = conn.find({"diff_reactome_enrich_id": ObjectId(reactome_enrich_id)}, target_cols)
    reactome_enrich_matrix = pd.DataFrame(list(reactome_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in reactome_enrich_matrix['seq_str']]
    reactome_enrich_matrix['seq_str'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    reactome_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export reactome_enrich matrix')
    return output

