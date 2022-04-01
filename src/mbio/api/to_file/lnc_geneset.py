# -*- coding: utf-8 -*-
# __author__ = 'sanger'

from __future__ import division
import os

import math

from biocluster.config import Config
from bson.objectid import ObjectId
import types
import json
import numpy as np
import re
from types import StringTypes
import gridfs
from collections import OrderedDict, defaultdict
import pandas as pd
from biocluster.file import getsize, exists
from biocluster.file import download
import shutil
from biocluster.api.file.lib.transfer import MultiFileTransfer
import sys

project_type = 'lnc_rna'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


# ---------------------基因集及相关gdq----------------------------------
def export_geneset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq 2019 04 修改
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    exp_id, geneset_id, geneset_type = data.split(";")
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
    seq_marker = {}
    geneset = set()
    if "all" not in geneset_id.lower():
        conn = db['sg_geneset_detail']
        for gs_id, gs_type in zip(re.split('\s*,\s*', geneset_id), re.split('\s*,\s*', geneset_type)):
            geneset_record = conn.find({"geneset_id": ObjectId(gs_id)})
            if not geneset_record:
                bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
            sub_geneset = geneset_record[0]['seq_list']
            geneset.update(sub_geneset)
            for seq_id in sub_geneset:
                seq_marker[seq_id] = gs_type

    # get all exp matrix
    conn = db['sg_exp_detail']
    if "refall" in geneset_id.lower():
        exp_records = conn.find({"exp_id": ObjectId(exp_id), 'is_new': False}, target_cols)
    else:
        exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: {}'.format(exp_id))
    exp_matrix = exp_matrix.set_index('seq_id')
    if geneset:
        needed_index = exp_matrix.index & geneset
        exp_matrix = exp_matrix.loc[needed_index, :]
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)

    marker = os.path.join(dir_path, 'seq_id_marker.json')
    with open(marker, 'w') as out_handler:
        json.dump(seq_marker, out_handler, indent=4)

    print('success to export expression matrix')
    return output


def export_geneset_from_query(data, option_name, dir_path, bind_obj=None):
    collection = db['sg_annotation_query_detail']
    main_collection = db['sg_annotation_query']
    my_result = main_collection.find_one({'task_id': data})
    if not my_result:
        bind_obj.set_error("意外错误，task_id:{}在sg_annotation_query中未找到！".format(data))
    if 'main_id' in my_result:
        query_id = my_result["main_id"]
    else:
        query_id = my_result['_id']
    results = collection.find({"query_id": ObjectId(query_id)})
    output = os.path.join(dir_path, "all_genes_transcripts.txt")
    with open(output, "w") as f:
        f.write('transcript_id\tgene_id\n')
        for result in results:
            transcript_id = ''
            if u'transcript_id' in result:
                transcript_id = result['transcript_id']
            gene_id = ''
            if u'gene_id' in result:
                gene_id = result['gene_id']
            f.write('{}\t{}\n'.format(transcript_id, gene_id))
    return output


def export_wgcna_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数仅为wgcna module分析服务。
    该函数会返回wgcna_prepare得到的exp_matrix路径，
    同时还会通过查询sg_annotation_query表导出geneid和genename的关系文件:seq_id2gene_name.txt
    """
    wgcna_prepare = db['sg_wgcna_prepare']
    prepare_id = data.strip()
    prepare_main = wgcna_prepare.find_one({"_id": ObjectId(prepare_id)})
    exp_matrix = prepare_main['exp_matrix']
    if re.match(r'^\w+://\S+/.+$', exp_matrix):
        download_s3_file(exp_matrix, os.path.join(bind_obj.work_dir, "exp_matrix.txt"))
        exp_matrix = os.path.join(bind_obj.work_dir, "exp_matrix.txt")
        if not os.path.exists(os.path.join(dir_path, "exp_matrix.txt")):
            os.link(exp_matrix, os.path.join(dir_path, "exp_matrix.txt"))
    else:
        if os.path.exists(os.path.join(dir_path, "exp_matrix.txt")):
            os.remove(os.path.join(dir_path, "exp_matrix.txt"))
        os.link(exp_matrix, os.path.join(dir_path, "exp_matrix.txt"))
    target_seqs = pd.read_table(exp_matrix, header=0, index_col=0).index
    task_id = prepare_main['task_id']
    annot_table = db['sg_annotation_query']
    try:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", "status": "end"})
    except:
        bind_obj.set_error("cannot find sg_annotation_query main table")
    else:
        if annot_main is None:
            annot_main = annot_table.find_one({"task_id": task_id, "type": "origin", "status": "end"})
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_annotation_query_detail']
    exp_level = bind_obj.sheet.option('exp_level')
    if exp_level[0].upper() == 'T':
        query_dict = dict(query_id=annot_main_id, )
        result_dict = dict(_id=0, gene_name=1, gene_id=1, transcript_id=1)
        result = annot_detail.find(query_dict, result_dict)
        gene2name = pd.DataFrame(list(result))
        gene2name.set_index('transcript_id', inplace=True)
    else:
        query_dict = dict(query_id=annot_main_id, is_gene=True)
        result_dict = dict(_id=0, gene_name=1, gene_id=1)
        result = annot_detail.find(query_dict, result_dict)
        gene2name = pd.DataFrame(list(result))
        gene2name.set_index('gene_id', inplace=True)
    gene2name = gene2name.loc[list(target_seqs), :]
    gene2name.reset_index(inplace=True)
    output = os.path.join(dir_path, "seq_id2gene_name.txt")
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    gene2name = pd.read_table(output, header=0)
    gene2name.fillna(method="pad", axis=1, inplace=True)
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    return exp_matrix + ';' + output


# added  by gdq for wgcna
def export_wgcna_relate_input(data, option_name, dir_path, bind_obj=None):
    module_id = data.strip()
    # export eigengenes
    eigengenes = db['sg_wgcna_module_eigengenes_detail']
    eigengenes_found = eigengenes.find({"module_id": ObjectId(module_id)}, {"_id": 0, "module_id": 0})
    eigengenes_pd = pd.DataFrame(list(eigengenes_found))
    eigengenes_pd.set_index("module", inplace=True)
    eigengenes_path = os.path.join(dir_path, "module_eigengenes.xls")
    eigengenes_pd.to_csv(eigengenes_path, sep='\t', header=True, index=True)
    # export exp matrix
    prepare_id = db['sg_wgcna_module'].find_one({"main_id": ObjectId(module_id)})['wgcna_prepare_id']
    exp_matrix = db['sg_wgcna_prepare'].find_one({"main_id": ObjectId(prepare_id)})['exp_matrix']
    # export each gene's module and gene_id/gene_name info
    exp_level = bind_obj.sheet.option('exp_level')
    membership = db['sg_wgcna_module_membership_detail']
    query_dict = dict(module_id=ObjectId(module_id))
    return_dict = dict(_id=0, seq_id=1, gene_name=1, module=1, kme=1, block_id=1)
    if exp_level == 'transcript':
        return_dict.update({"gene_id": 1})
    result = membership.find(query_dict, return_dict)
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("seq_id", inplace=True)
    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    return exp_matrix + ';' + eigengenes_path + ";" + gene_annot


def get_after_qc_bam_path(data, option_name, dir_path, bind_obj=None):
    db = Config().get_mongo_client(mtype="lnc_rna")[Config().get_mongo_dbname("lnc_rna")]
    conn = db['sg_specimen']
    result = conn.find({"task_id": data, "about_qc": "after"})
    sorted_result = sorted(result, key=lambda k: k['_id'])
    output = os.path.join(dir_path, option_name)
    with open(output, 'w') as f:
        for i in sorted_result:
            f.write(i['bam_path'] + "\n")
    return output


# ---------------------从蛋白复制过来的------------------------

def export_compare_exp_fc(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析accession_id 和 fc
    '''

    diff_id = bind_obj.sheet.option('diff_id')
    compare_group = bind_obj.sheet.option('compare_group')
    target_cols = OrderedDict(seq_id=1, log2fc=1, _id=0)

    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id), "compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    output = os.path.join(dir_path, option_name)
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export expression matrix')
    return output


def export_go_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    go_enrich_id = bind_obj.sheet.option('enrich_id')
    go_type = bind_obj.sheet.option('go_type')

    target_cols = OrderedDict(go_id=1, go_type=1, discription=1, p_corrected=1, p_uncorrected=1, seq_list=1, depth=1,
                              _id=0)
    bind_obj.logger.debug("导出GO富集表")
    # get geneset
    conn = db['sg_geneset_go_enrich_detail']
    if go_type == 'ALL' or go_type == 'All':
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id)}, target_cols)
    else:
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id), "go_type": go_type}, target_cols)
    go_enrich_matrix = pd.DataFrame(list(go_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in go_enrich_matrix['seq_list']]
    go_enrich_matrix['seq_list'] = seq_list2
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
    target_cols = OrderedDict(_id=0, id=1, term=1, pvalue=1, corrected_pvalue=1, seq_list=1, kegg_type=1)
    bind_obj.logger.debug("导出KEGG富集表")
    # get geneset
    conn = db['sg_geneset_kegg_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(kegg_enrich_id))
    kegg_enrich_records = conn.find({"kegg_enrich_id": ObjectId(kegg_enrich_id)}, target_cols)
    kegg_enrich_matrix = pd.DataFrame(list(kegg_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in kegg_enrich_matrix['seq_list']]
    kegg_enrich_matrix['seq_list'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    kegg_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export kegg_enrich matrix')
    return output


def export_gene_list_ppi(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "%s.txt" % option_name)
    bind_obj.logger.debug("正在导出蛋白集")
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(data)))
    results = collection.find_one({"geneset_id": ObjectId(data)})["seq_list"]
    with open(gene_list_path, "wb") as f:
        f.write("accession_id" + "\n")
        for result in results:
            f.write(result + "\n")
    bind_obj.logger.debug("蛋白集导出成功！")
    return gene_list_path


def download_s3_file(path, to_path):
    """
    判断文件是否在对象存储上
    """
    if not to_path.startswith("/"):
        to_path = os.path.join(self.work_dir, to_path)
    if os.path.exists(to_path):
        os.remove(to_path)
    elif os.path.exists(path):
        to_path = path
    elif exists(path):
        download(path, to_path)
    else:
        print
        'file can not find'
    return to_path


# added for transcription factor analysis  --------------------------------------------
def get_all_pep_seq(data, option_name, dir_path, bind_obj=None):
    pep_db_path = data.strip()
    result_path = os.path.join(dir_path, "all_pep.fa")
    if re.match(r'^\w+://\S+/.+$', pep_db_path) or re.match(r'/mnt/ilustre', pep_db_path):
        transfer = MultiFileTransfer()
        transfer.add_download(pep_db_path, bind_obj.work_dir + "/")
        transfer.perform()
        # download_s3_file(pep_db_path, os.path.join(bind_obj.work_dir, "pep.db"))
        pep_db_path = os.path.join(bind_obj.work_dir, "refrna_seqs.db")
    bind_obj.logger.debug("文件为 {}".format(pep_db_path))
    import sqlite3
    conn = sqlite3.connect(pep_db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT transcript_id,pep_seq FROM trans_annot")
    with open(result_path, 'w') as fw:
        # pep_id 就是转录本的id
        for pep_id, pep_seq in cursor.fetchall():
            fw.write('>{}\n{}\n'.format(pep_id, pep_seq))
    return result_path


def get_gene_detail(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip()
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", 'status': 'end'})
    if not annot_main:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    if not annot_main:
        task_table = db['sg_task']
        task_main = task_table.find_one({"task_id": task_id})
        if os.path.exists(task_main['annot']):
            gene_annot = os.path.join(dir_path, "seq_annot.xls")
            os.system("cut -f 1-4 " + task_main['annot'] + " |awk '{print $2\"\\t\"$1\"\\t\"$4}' > " + gene_annot)
            return gene_annot
            # os.system("cut -f 1-4 {} |awk '{print $2\"\t\"$1\"\t\"$4}' > {} ".format(task_main['annot'], gene_annot))
        else:
            bind_obj.set_error("Not Found in sg_annotation_query by query {}".format(task_id))
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_annotation_query_detail']
    query_dict = dict(query_id=annot_main_id, )
    result_dict = dict(_id=0, gene_name=1, gene_id=1, description=1, transcript_id=1)
    result = annot_detail.find(query_dict, result_dict)
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("transcript_id", inplace=True)
    result_pd = result_pd.loc[:, ["gene_id", "gene_name", "description"]]
    result_pd.columns = ["gene_id", "gene_name", "gene_desc"]
    # result_pd = result_pd.loc[:, ["gene_id", "gene_name"]]
    # result_pd.columns = ["gene_id", "gene_name"]
    result_pd["gene_name"][result_pd['gene_name'] == "-"] = None
    result_pd = result_pd.fillna(method="pad", axis=1)
    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    os.system(
        r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot


# ---------------------基因集及相关zhaozhipeng----------------------------------
def export_genesets_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    s_list = data.split(";")
    exp_id = s_list[0]
    geneset_ids = s_list[1:]
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
    geneset = set()
    bind_obj.logger.debug('=========== %s ++++++++++++++' % ','.join(geneset_ids))
    for geneset_id in geneset_ids:
        if "all" not in geneset_id.lower():
            conn = db['sg_geneset_detail']
            geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
            if not geneset_records:
                bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
            geneset.update(geneset_records['seq_list'])

    # get all exp matrix
    conn = db['sg_exp_detail']
    exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: {}'.format(exp_id))

    exp_matrix = exp_matrix.set_index('seq_id')
    if geneset:
        indexs = exp_matrix.index & {i for i in geneset}
        exp_matrix = exp_matrix.loc[indexs, :]
    exp_matrix.dropna(how='all')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output


def export_genes_list(data, option_name, dir_path, bind_obj=None):
    """

    :param data:
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    bind_obj.logger.debug('=========== %s ++++++++++++++' % data)
    geneset_ids = [item.split(',') for item in data.split(";")]
    geneset = {}
    for gene_type, geneset_id in geneset_ids:
        if "all" not in geneset_id.lower():
            conn = db['sg_geneset_detail']
            geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
            if not geneset_records:
                bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
            for i in geneset_records['seq_list']:
                geneset[i] = 'lncRNA' if gene_type == 'l_geneset_id' else 'mRNA'

    output = os.path.join(dir_path, 'genes_info.json')
    with open(output, 'w') as in_handler:
        json.dump(geneset, in_handler, indent=4)

    return output


def export_diff_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """

    :param data: {'diff_group': 'xxx|yyy', 'sg_diff': '5c9da1ba17b2bf206a55d350'}
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    bind_obj.logger.debug(' :: export diff_exp.txt file')
    compare = data['diff_group']
    sg_diff_id = data['sg_diff_id']

    conn = db['sg_diff_detail']
    recoreds = conn.find({'diff_id': ObjectId(sg_diff_id), 'compare': compare})
    res_data = []
    for dic in recoreds:
        tmp_dic = {'seq_id': dic['seq_id'], 'log2fc': dic['log2fc'], 'regulate': dic['regulate'],
                   'avg_exp': (dic['group1'] + dic['group2']) / 2, 'significant': dic['significant']}
        res_data.append(tmp_dic)
    result_df = pd.DataFrame(res_data)
    output = os.path.join(dir_path, 'diff_exp.txt')
    result_df.to_csv(output, sep='\t', index=False)
    bind_obj.logger.debug(' :: export diff_exp.txt file :: completed')

    return output


def export_rna_type_matrix(data, option_name, dir_path, bind_obj=None):
    """

    :param data: {'task_id': 'lnc_rna', 'exp_level': 'T', 'sg_diff': '5c9da1ba17b2bf206a55d350'}
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    bind_obj.logger.debug(' :: export seq_id_type.txt file')
    task_id = data['task_id']
    sg_diff_id = data['diff_exp']
    exp_level = data['exp_level']
    # {'task_id': data.task_id, 'exp_level': 'T', 'diff_exp': data.diff_exp}
    conn = db['sg_diff']
    bind_obj.logger.debug(str({'main_id': ObjectId(sg_diff_id), 'exp_level': exp_level, 'task_id': task_id}))
    # recored = conn.find_one({'main_id': ObjectId(sg_diff_id), 'exp_level': exp_level, 'task_id': task_id})
    # 因为main_id本身唯一,所以删除多余查找字段
    recored = conn.find_one({'main_id': ObjectId(sg_diff_id)})
    exp_id = ObjectId(recored['exp_id'])

    conn = db['sg_exp_detail']
    recoreds = conn.find({'exp_id': exp_id}, {'seq_id': 1, 'rna_type': 1, '_id': 0})
    result_df = pd.DataFrame([i for i in recoreds])
    output = os.path.join(dir_path, 'seq_id_type.txt')
    result_df.to_csv(output, sep='\t', index=False)
    bind_obj.logger.debug(' :: export seq_id_type.txt file :: completed')

    return output


def export_cis_matrix2(data, option_name, dir_path, bind_obj=None):
    """

    :param data: {'task_id': 'lnc_rna', 'exp_level': 'T', 'sg_diff': '5c9da1ba17b2bf206a55d350'}
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    bind_obj.logger.debug(' :: export cis_targets.txt file')

    # 获取基因转录本对应关系， 只要一组

    task_id = bind_obj.sheet.option("task_id")
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "origin", "status": "end"})
    annot_detail = db['sg_annotation_query_detail']
    result = annot_detail.find({"query_id": annot_main["main_id"], "is_gene": True},
                               {"_id": 0, "transcript_id": 1, "gene_id": 1})
    g2t_df = pd.DataFrame(list(result))

    conn = db['sg_target_cistrans_detail']
    recoreds = conn.find({'target_cistrans_id': ObjectId(data), "target_type": "cis"},
                         {'lncrna_id': 1, 'gene_id': 1, 'mgene_id': 1, '_id': 0})
    data = [i for i in recoreds]
    output = os.path.join(dir_path, 'cis_targets.txt')
    if data:
        result_df = pd.DataFrame(data)
        result_df = result_df.merge(g2t_df, left_on = 'mgene_id', right_on = 'gene_id', suffixes=("", "_ri"))
        result_df.drop(columns=["gene_id_ri"],  inplace=True)
        result_df.rename(columns={'transcript_id': "mrna_id"}, inplace=True)
        result_df.to_csv(output, sep='\t', index=False)
    else:
        bind_obj.logger.debug('cis_targets.txt 文件为空, %s' % len(data))
        with open(output, 'w') as out_handler:
            out_handler.write('lncrna_id\tgene_id\tmrna_id\tmgene_id\n')
    bind_obj.logger.debug(' :: export cis_targets.txt file :: completed')

    return output

def export_cis_matrix(data, option_name, dir_path, bind_obj=None):
    """

    :param data: {'task_id': 'lnc_rna', 'exp_level': 'T', 'sg_diff': '5c9da1ba17b2bf206a55d350'}
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    bind_obj.logger.debug(' :: export cis_targets.txt file')

    conn = db['sg_target_cis_detail']
    recoreds = conn.find({'target_cis_id': ObjectId(data)},
                         {'lncrna_id': 1, 'gene_id': 1, 'mrna_id': 1, 'mgene_id': 1, '_id': 0})
    data = [i for i in recoreds]
    output = os.path.join(dir_path, 'cis_targets.txt')
    if data:
        result_df = pd.DataFrame(data)

        result_df.to_csv(output, sep='\t', index=False)
    else:
        bind_obj.logger.debug('cis_targets.txt 文件为空, %s' % len(data))
        with open(output, 'w') as out_handler:
            out_handler.write('lncrna_id\tgene_id\tmrna_id\tmgene_id\n')
    bind_obj.logger.debug(' :: export cis_targets.txt file :: completed')

    return output

def export_trans_matrix2(data, option_name, dir_path, bind_obj=None):
    """

    :param data: {'task_id': 'lnc_rna', 'exp_level': 'T', 'sg_diff': '5c9da1ba17b2bf206a55d350'}
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    bind_obj.logger.debug(' :: export trans_targets.txt file')

    # 获取基因转录本对应关系， 只要一组

    task_id = bind_obj.sheet.option("task_id")
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "origin", "status": "end"})
    annot_detail = db['sg_annotation_query_detail']
    result = annot_detail.find({"query_id": annot_main["main_id"], "is_gene": True},
                               {"_id": 0, "transcript_id": 1, "gene_id": 1})
    g2t_df = pd.DataFrame(list(result))

    conn = db['sg_target_cistrans_detail']
    recoreds = conn.find({'target_cistrans_id': ObjectId(data), "target_type": "trans"},
                         {'lncrna_id': 1, 'gene_id': 1, 'mgene_id': 1, '_id': 0})

    data = []
    check = defaultdict(int)
    for recored in recoreds:
        lnc_id = recored['lncrna_id']
        num = check[lnc_id]
        if num >= 3:
            continue
        check[lnc_id] += 1
        data.append(recored)

    output = os.path.join(dir_path, 'trans_targets.txt')
    if data:
        result_df = pd.DataFrame(data)
        result_df = result_df.merge(g2t_df, left_on = 'mgene_id', right_on = 'gene_id', suffixes=("", "_ri"))
        result_df.drop(columns=["gene_id_ri"],  inplace=True)
        result_df.rename(columns={'transcript_id': "mrna_id"}, inplace=True)

        result_df.to_csv(output, sep='\t', index=False)
    else:
        bind_obj.logger.debug('trans_targets.txt 文件为空, %s' % len(data))
        with open(output, 'w') as out_handler:
            out_handler.write('lncrna_id\tgene_id\tmrna_id\tmgene_id\n')
    bind_obj.logger.debug(' :: export trans_targets.txt file :: completed')

    return output



def export_trans_matrix(data, option_name, dir_path, bind_obj=None):
    """

    :param data: {'task_id': 'lnc_rna', 'exp_level': 'T', 'sg_diff': '5c9da1ba17b2bf206a55d350'}
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    bind_obj.logger.debug(' :: export trans_targets.txt file')

    conn = db['sg_target_trans_detail']
    recoreds = conn.find({'target_trans_id': ObjectId(data)},
                         {'lncrna_id': 1, 'gene_id': 1, 'mrna_id': 1, 'mgene_id': 1, '_id': 0})
    data = []
    check = defaultdict(int)
    for recored in recoreds:
        lnc_id = recored['lncrna_id']
        num = check[lnc_id]
        if num >= 3:
            continue
        check[lnc_id] += 1
        data.append(recored)

    output = os.path.join(dir_path, 'trans_targets.txt')
    if data:
        result_df = pd.DataFrame(data)
        result_df.to_csv(output, sep='\t', index=False)
    else:
        bind_obj.logger.debug('trans_targets.txt 文件为空, %s' % len(data))
        with open(output, 'w') as out_handler:
            out_handler.write('lncrna_id\tgene_id\tmrna_id\tmgene_id\n')
    bind_obj.logger.debug(' :: export trans_targets.txt file :: completed')

    return output


def export_multi_gene_list(data, option_name, dir_path, bind_obj=None):
    data = json.loads(data)
    geneset_id_list = data['geneset_id'].split(",")
    source = data['source']
    multi_geneset_path = dir_path + "/multi_geneset_list"
    main_collection = db['sg_geneset']
    collection = db['sg_geneset_detail']
    with open(multi_geneset_path, "wb") as out_handler:
        if len(geneset_id_list) == 1 and source == 'diff_exp':
            geneset_id = geneset_id_list[0]
            my_result = main_collection.find_one({'main_id': ObjectId(geneset_id)})
            if not my_result:
                bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(geneset_id)))
            geneset_name = my_result["name"]
            results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
            seq_list = results["seq_list"]
            if 'regulate_list' not in results:
                raise bind_obj.set_error(
                    'geneset_detail table must contain "regulate_list" field when provide one geneset id')
            regulate_list = results["regulate_list"]
            up_list, down_list = [], []
            for seq_id, regulate in zip(seq_list, regulate_list):
                if regulate == 'up':
                    up_list.append(seq_id)
                elif regulate == 'down':
                    down_list.append(seq_id)
                else:
                    continue
            out_handler.write(geneset_name + '_up\t' + ','.join(up_list) + '\n')
            out_handler.write(geneset_name + '_down\t' + ','.join(down_list) + '\n')
        else:
            for n, gi in enumerate(geneset_id_list):
                my_result = main_collection.find_one({'main_id': ObjectId(gi)})
                if not my_result:
                    bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(gi)))
                out_handler.write(my_result["name"] + "\t")
                results = collection.find_one({"geneset_id": ObjectId(gi)})
                out_handler.write(",".join(results["seq_list"]) + "\n")

    return multi_geneset_path


def export_go_detail(data, option_name, dir_path, bind_obj=None):
    geneset_id_list = data.split(",")
    multi_geneset_path = dir_path + "/multi_geneset_list"
    main_collection = db['sg_geneset']
    collection = db['sg_geneset_detail']
    with open(multi_geneset_path, "wb") as out_handler:
        if len(geneset_id_list) == 2:
            for n, gi in enumerate(geneset_id_list):
                my_result = main_collection.find_one({'main_id': ObjectId(gi)})
                if not my_result:
                    bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(gi)))
                out_handler.write(my_result["name"] + "\t")
                results = collection.find_one({"geneset_id": ObjectId(gi)})
                out_handler.write(",".join(results["seq_list"]) + "\n")
        else:
            geneset_id = geneset_id_list[0]
            my_result = main_collection.find_one({'main_id': ObjectId(geneset_id)})
            if not my_result:
                bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(geneset_id)))
            geneset_name = my_result["name"]
            results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
            seq_list = results["seq_list"]
            if 'regulate_list' not in results:
                raise bind_obj.set_error(
                    'geneset_detail table must contain "regulate_list" field when provide one geneset id')
            regulate_list = results["regulate_list"]
            up_list, down_list = [], []
            for seq_id, regulate in zip(seq_list, regulate_list):
                if regulate == 'up':
                    up_list.append(seq_id)
                else:
                    down_list.append(seq_id)
            out_handler.write(geneset_name + '_up\t' + ','.join(up_list) + '\n')
            out_handler.write(geneset_name + '_down\t' + ','.join(down_list) + '\n')

    return multi_geneset_path


def export_kegg_detail(data, option_name, dir_path, bind_obj=None):
    """

    :param data: {''}
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    geneset_id_list = data.split(",")
    multi_geneset_path = dir_path + "/multi_geneset_list"
    main_collection = db['sg_geneset']
    collection = db['sg_geneset_detail']
    with open(multi_geneset_path, "wb") as out_handler:
        if len(geneset_id_list) == 2:
            for n, gi in enumerate(geneset_id_list):
                my_result = main_collection.find_one({'main_id': ObjectId(gi)})
                if not my_result:
                    bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(gi)))
                out_handler.write(my_result["name"] + "\t")
                results = collection.find_one({"geneset_id": ObjectId(gi)})
                out_handler.write(",".join(results["seq_list"]) + "\n")
        else:
            geneset_id = geneset_id_list[0]
            my_result = main_collection.find_one({'main_id': ObjectId(geneset_id)})
            if not my_result:
                bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(geneset_id)))
            geneset_name = my_result["name"]
            results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
            seq_list = results["seq_list"]
            if 'regulate_list' not in results:
                bind_obj.set_error(
                    'geneset_detail table must contain "regulate_list" field when provide one geneset id')
            regulate_list = results["regulate_list"]
            up_list, down_list = [], []
            for seq_id, regulate in zip(seq_list, regulate_list):
                if regulate == 'up':
                    up_list.append(seq_id)
                else:
                    down_list.append(seq_id)
            out_handler.write(geneset_name + '_up\t' + ','.join(up_list) + '\n')
            out_handler.write(geneset_name + '_down\t' + ','.join(down_list) + '\n')

    return multi_geneset_path


def export_enrich_heatmap(data, option_name, dir_path, bind_obj=None):
    """
        ['task_id', 'submit_location', 'geneset_type', 'enrich_type',
        'geneset_enrich_ids', 'go_level', 'task_type', 'pvalue_type', 'top_num', 'ids_list', 'go_level']

    :param data:
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    ids_filer_list = None
    if 'ids_list' in data:
        ids_filer_list = {i for i in data['ids_list']}
    query_dict = {}
    filter_dict = dict(_id=0)
    out_file = os.path.join(dir_path, 'heatmap_data.xls')

    if data['enrich_type'].lower() == 'kegg':
        if data['pvalue_type'].lower() == 'pvalue':
            # if ids_filer_list is None:
            #     query_dict['pvalue'] = {'$lte': data['p_value']}
            filter_dict['pvalue'] = 1
            p_field = 'pvalue'
        else:
            # if ids_filer_list is None:
            #     query_dict['corrected_pvalue'] = {'$lte': data['p_value']}
            filter_dict['corrected_pvalue'] = 1
            p_field = 'corrected_pvalue'

        main_collection = db['sg_geneset_kegg_enrich']
        collection = db['sg_geneset_kegg_enrich_detail']
        ids = [ObjectId(i) for i in data['geneset_enrich_ids'].strip().split(',')]
        filter_dict['id'] = 1
        filter_dict['kegg_type'] = 1
        filter_dict['discription'] = 1

        data_list = {}
        detail_info = defaultdict(dict)
        needed_ids = set()
        geneset_names= []
        min_non_zero = 1
        has_zero = False
        for id_obj in ids:
            query_dict['kegg_enrich_id'] = id_obj
            name = main_collection.find_one({'_id': id_obj})['geneset_name']
            geneset_names.append(name)
            sub_data = {}
            data_list[name] = sub_data
            for dic in collection.find(query_dict, filter_dict):
                kegg_id = dic['id']
                sub_dic = detail_info[kegg_id]
                if 'kegg_type' not in sub_dic:
                    sub_dic['kegg_type'] = dic.pop('kegg_type')
                    sub_dic['discription'] = dic.pop('discription')
                p_value = dic[p_field]
                sub_data[kegg_id] = p_value
                if p_value == 0:
                    has_zero = True
                else:
                    min_non_zero = min(p_value, min_non_zero)
            if ids_filer_list is None:
                # 过滤top num 的ids
                top_num_ids = {t for t in
                               sorted(((k, v) for k, v in sub_data.items()), key=lambda items: items[1])[
                               : data['top_num']]}
                # 过滤p_value值
                p_filter_ids = {t[0] for t in top_num_ids if t[1] <= data['p_value']}
                # 取所有ids的并集
                needed_ids.update(p_filter_ids)

        df = pd.DataFrame(data_list)
        ids_filer_list = ids_filer_list or needed_ids
        if len(ids_filer_list) == 0:
            bind_obj.set_error('提供的id集合和基因集没有交集，请重新输入参数')
        if len(geneset_names) != len(set(geneset_names)):
            bind_obj.set_error('提供的id集合名称重复：' + str(geneset_names))

        df = df.loc[df.index & ids_filer_list]
        if has_zero:  # 0 值取最小正实数的十分之一 [比最小正实数小一个数量级]
            df[df == 0] = min_non_zero / 10
        df.fillna(1, inplace=True)
        df = -np.log10(df)
        df[df == 0] = 0
        indexs = df.index
        df['kegg_type'] = [detail_info[i]['kegg_type'] for i in indexs]
        df['discription'] = [detail_info[i]['discription'] for i in indexs]

        df.index.name = 'kegg_id'
        df.to_csv(out_file, sep='\t', index=True, header=True)
    else:
        if data['pvalue_type'].lower() == 'pvalue':
            # if ids_filer_list is None:
            #     query_dict['p_uncorrected'] = {'$lte': data['p_value']}
            filter_dict['p_uncorrected'] = 1
            p_field = 'p_uncorrected'
        else:
            # if ids_filer_list is None:
            #     query_dict['p_corrected'] = {'$lte': data['p_value']}
            filter_dict['p_corrected'] = 1
            p_field = 'p_corrected'

        main_collection = db['sg_geneset_go_enrich']
        collection = db['sg_geneset_go_enrich_detail']
        ids = [ObjectId(i) for i in data['geneset_enrich_ids'].strip().split(',')]
        filter_dict['go_id'] = 1
        filter_dict['go_type'] = 1
        filter_dict['discription'] = 1
        if 'go_level' in data:
            go_level = data['go_level']
            if go_level.lower() != 'all':
                query_dict['go_type'] = go_level.upper()

        data_list = {}
        detail_info = defaultdict(dict)
        needed_ids = set()
        geneset_names= []
        min_non_zero = 1
        has_zero = False
        for id_obj in ids:
            query_dict['go_enrich_id'] = id_obj
            name = main_collection.find_one({'_id': id_obj})['geneset_name']
            geneset_names.append(name)
            sub_data = {}
            data_list[name] = sub_data
            print(query_dict, '========================')
            for dic in collection.find(query_dict, filter_dict):
                go_id = dic['go_id']
                sub_dic = detail_info[go_id]
                if 'go_type' not in sub_dic:
                    sub_dic['go_type'] = dic.pop('go_type')
                    sub_dic['discription'] = dic.pop('discription')
                p_value = dic[p_field]
                sub_data[go_id] = p_value
                if p_value == 0:
                    has_zero = True
                else:
                    min_non_zero = min(p_value, min_non_zero)

            if ids_filer_list is None:
                # 过滤top num 的ids
                top_num_ids = {t for t in
                               sorted(((k, v) for k, v in sub_data.items()), key=lambda items: items[1])[
                               : data['top_num']]}
                # 过滤p_value值
                p_filter_ids = {t[0] for t in top_num_ids if t[1] <= data['p_value']}
                # 取所有ids的并集
                needed_ids.update(p_filter_ids)

        df = pd.DataFrame(data_list)
        ids_filer_list = ids_filer_list or needed_ids
        if len(ids_filer_list) == 0:
            bind_obj.set_error('提供的id集合和基因集没有交集，请重新输入参数')
        df = df.loc[df.index & ids_filer_list]
        if has_zero:  # 0 值取最小正实数的十分之一
            df[df == 0] = min_non_zero / 10
        df.fillna(1, inplace=True)
        df = -np.log10(df)
        df[df == 0] = 0
        indexs = df.index
        df['go_type'] = [detail_info[i]['go_type'] for i in indexs]
        df['discription'] = [detail_info[i]['discription'] for i in indexs]
        df.fillna(1, inplace=True)

        df.index.name = 'go_id'
        df.to_csv(out_file, sep='\t', index=True, header=True)

    return out_file
