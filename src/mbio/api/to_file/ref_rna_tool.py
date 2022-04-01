# -*- coding: utf-8 -*-
# __author__ = 'sanger'

from biocluster.config import Config
import os
import json
from bson.objectid import ObjectId
from collections import OrderedDict
import sys
import pandas as pd


project_type = 'ref_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

def checkwargs(**kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('start checking to_file arguments')
        for k, v in kwargs.items():
            kwargs['bind_obj'].logger.debug('{} = {}'.format(k, v))

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

def export_circularbar(data, option_name, dir_path, bind_obj=None):
    if data.split('____')[2] == "geneset_go":
        bind_obj.logger.debug("正在导出GO分类信息")
        collection = db['sg_geneset_go_class_detail']
        main_collection = db['sg_geneset_go_class']
        go_id = data.split('____')[0]
        geneset = data.split('____')[1] + "_num"
        my_result = main_collection.find_one({'main_id': ObjectId(go_id)})
        if not my_result:
            bind_obj.set_error("意外错误，main_id:{}在sg_geneset_go_class中未找到！".format(ObjectId(go_id)))
        target_cols = OrderedDict(zip(['go','go_type',geneset,'_id'], [1,1,1,0]))
        results = collection.find({"go_regulate_id": ObjectId(go_id)}, target_cols)
        go_df = pd.DataFrame(list(results))
        columnc_order = ['go', 'go_type', geneset]
        go_df = go_df[columnc_order]
        # sort the matrix by seq_num and goterm_2 
        go_df[[geneset]] = go_df[[geneset]].astype(int)
        go_df = go_df.sort_values(by=['go_type',geneset], ascending=[True,False])
        go_df = pd.concat([i[1][:20] for i in go_df.groupby(by=["go_type"])])
        go_df = go_df.rename(columns={'go': 'Go', 'go_type': 'Go Type',geneset: 'Number'})
        go_df = go_df[go_df['Number']>0]
        #
        output = os.path.join(dir_path, "circularbar.txt")
        go_df.to_csv(output, sep='\t', header=True, index=False)
        print('success to export go dataframe matrix')
        return output
    elif data.split('____')[2] == "geneset_kegg":
        bind_obj.logger.debug("正在导出KEGG分类信息")
        collection = db['sg_geneset_kegg_class_detail']
        main_collection = db['sg_geneset_kegg_class']
        kegg_id = data.split('____')[0]
        geneset = data.split('____')[1] + "_numbers"
        my_result = main_collection.find_one({'main_id': ObjectId(kegg_id)})
        if not my_result:
            bind_obj.set_error("意外错误，main_id:{}在sg_geneset_kegg_class中未找到！".format(ObjectId(kegg_id)))
        target_cols = OrderedDict(zip(['second_category','first_category',geneset,'_id'], [1,1,1,0]))
        results = collection.find({"kegg_id": ObjectId(kegg_id)}, target_cols)
        kegg_df = pd.DataFrame(list(results))
        columnc_order = ['second_category','first_category', geneset]
        kegg_df = kegg_df[columnc_order]
        # sort the matrix by seq_num and goterm_2 
        kegg_df[[geneset]] = kegg_df[[geneset]].astype(int)
        kegg_df = kegg_df.sort_values(by=['first_category',geneset], ascending=[True,False])
        # kegg_df = pd.concat([i[1][:20] for i in kegg_df.groupby(by=["kegg_type"])])
        kegg_df = kegg_df.rename(columns={'second_category': 'Second Category', 'first_category': 'First Category',geneset: 'Number'})
        kegg_df = kegg_df[kegg_df['Number']>0]
        kegg_df = kegg_df.drop_duplicates(subset=['Second Category'],keep='first')
        #
        output = os.path.join(dir_path, "circularbar.txt")
        kegg_df.to_csv(output, sep='\t', header=True, index=False)
        print('success to export go dataframe matrix')
        return output
    else:
        bind_obj.set_error("意外错误，对{}进行tofie时未发现可判断的web_page".format(data))
