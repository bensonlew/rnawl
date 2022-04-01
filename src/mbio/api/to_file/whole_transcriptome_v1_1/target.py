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



def export_geneset_lncset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    if len(data.split(";")) == 5:
        exp_id, geneset_id, lncset_id, level, is_rmbe = data.split(";")
        print "is is {} | {} | {}".format(exp_id, geneset_id, lncset_id)
    else:
        print "gene set wrong"
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    # export group info
    with open(dir_path+'/group_info.txt', 'w') as f:
        f.write('#sample\tgroup\n')
        group_id = bind_obj.sheet.option('group_id').lower()
        for g in group_dict:
            for s in group_dict[g]:
                if group_id == "all":
                    g = s
                f.write('{}\t{}\n'.format(s, g))
    # decide output columns
    if level == "G":
        target_cols = OrderedDict(gene_id=1, _id=0)
    else:
        target_cols = OrderedDict(transcript_id=1, _id=0)
    if is_rmbe == 'false':
        for each in samples:
            target_cols[each] = 1
    if is_rmbe == 'true':
        for each in samples:
            target_cols['{}_batch'.format(each)] = 1
    # get geneset
    if geneset_id == "":
        geneset = list()
    elif "all" not in geneset_id.lower():
        conn = db['geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        if not geneset_records:
            bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
        geneset = geneset_records['seq_list']
    else:
        geneset = list()

    # get_lncrnaset
    if lncset_id == "":
        lncset = list()
    elif "all" not in lncset_id.lower():
        conn = db['geneset_detail']
        lncset_records = conn.find_one({"geneset_id": ObjectId(lncset_id)})
        if not lncset_records:
            bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
        lncset = lncset_records['seq_list']
    else:
        lncset = list()
    # get all exp matrix
    conn = db['exp_detail']
    if "refall" in geneset_id.lower():
        exp_records = conn.find({"exp_id": ObjectId(exp_id), 'kind': "ref", 'category': "mRNA"}, target_cols)
    else:
        exp_records = conn.find({"exp_id": ObjectId(exp_id), 'category': "mRNA"}, target_cols)

    if "refall" in lncset_id.lower():
        lnc_exp_records = conn.find({"exp_id": ObjectId(exp_id), 'kind': "ref", 'category': "lncRNA"}, target_cols)
    else:
        lnc_exp_records = conn.find({"exp_id": ObjectId(exp_id), 'category': "lncRNA"}, target_cols)

    exp_matrix = pd.DataFrame(list(exp_records))
    '''
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: {}'.format(exp_id))
    '''

    lnc_exp_matrix = pd.DataFrame(list(lnc_exp_records))
    '''
    if lnc_exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: {}'.format(exp_id))
    '''


    if level == "G":
        exp_matrix.to_csv(dir_path + "/test1.tsv", sep="\t")
        exp_matrix.rename(columns={"gene_id": "seq_id"}, inplace = True)
        if is_rmbe == 'true':
            exp_matrix.rename(columns=dict(zip(['{}_batch'.format(i) for i in samples], samples)), inplace=True)
        if is_rmbe == 'false':
            pass
        lnc_exp_matrix.rename(columns={"gene_id": "seq_id"}, inplace = True)
    else:
        exp_matrix.to_csv(dir_path + "/test2.tsv", sep="\t")
        exp_matrix.rename(columns = {"transcript_id": "seq_id"}, inplace = True)
        lnc_exp_matrix.rename(columns = {"transcript_id": "seq_id"}, inplace = True)




    print "geneset is {}".format(geneset)
    print "lncset is {}".format(lncset)

    if len(lncset) != 0  or "all" in lncset_id.lower():
        lnc_exp_matrix = lnc_exp_matrix.set_index('seq_id')
        if len(lncset) != 0:
            lnc_exp_matrix = lnc_exp_matrix.loc[lncset, :]

    if len(geneset) != 0  or "all" in geneset_id.lower() :
        exp_matrix = exp_matrix.set_index('seq_id')
        if len(geneset) != 0:
            exp_matrix = exp_matrix.loc[geneset, :]

    if geneset_id != "" :
        exp_matrix = exp_matrix
        if lncset_id != "" :
            exp_matrix = exp_matrix.append(lnc_exp_matrix)
    else:
        exp_matrix = lnc_exp_matrix



    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_geneset_from_query(data, option_name, dir_path, bind_obj=None):
    collection = db['annotation_query_detail']
    main_collection = db['annotation_query']
    my_result = main_collection.find_one({'task_id': data})
    if not my_result:
        bind_obj.set_error("意外错误，task_id:{}在annotation_query中未找到！".format(data))
    if 'main_id' in my_result:
        query_id = my_result["main_id"]
    else:
        query_id = my_result['_id']
    results = collection.find({"query_id": ObjectId(query_id)})
    output = os.path.join(dir_path, "all_genes_transcripts.txt")
    with open(output, "w") as f:
        f.write('transcript_id\tgene_id\n')
        for result in results:
            transcript_id=''
            if u'transcript_id' in result:
                transcript_id=result['transcript_id']
            gene_id=''
            if u'gene_id' in result:
                gene_id = result['gene_id']
            f.write('{}\t{}\n'.format(transcript_id, gene_id))
    return output


def export_wgcna_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数仅为wgcna module分析服务。
    该函数会返回wgcna_prepare得到的exp_matrix路径，
    同时还会通过查询annotation_query表导出geneid和genename的关系文件:seq_id2gene_name.txt
    """
    wgcna_prepare = db['wgcna_prepare']
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
    annot_table = db['annotation_query']
    lnc_annot_table = db['known_lncrna_identify']
    lncnew_annot_table = db['new_lncrna_predict']
    try:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        lnc_annot_main = lnc_annot_table.find_one({"task_id": task_id, "status": "end"})
        lncnew_annot_main = lncnew_annot_table.find_one({"task_id": task_id, "status": "end"})
    except:
        bind_obj.set_error("cannot find annotation_query main table")
    else:
        if annot_main is None:
            annot_main = annot_table.find_one({"task_id": task_id, "type": "origin", "status": "end"})
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']

    if "main_id" not in lnc_annot_main:
        lnc_annot_main_id = lnc_annot_main['_id']
    else:
        lnc_annot_main_id = lnc_annot_main['main_id']

    if "main_id" not in lncnew_annot_main:
        lncnew_annot_main_id = lncnew_annot_main['_id']
    else:
        lncnew_annot_main_id = lncnew_annot_main['main_id']

    annot_detail = db['annotation_query_detail']
    lnc_annot_detail = db['known_lncrna_identify_detail']
    lncnew_annot_detail = db['new_lncrna_predict_detail']

    exp_level = bind_obj.sheet.option('exp_level')
    if exp_level[0].upper() == 'T':
        query_dict = dict(query_id=annot_main_id,)
        result_dict = dict( _id=0, gene_name=1, gene_id=1, transcript_id=1)
        result = annot_detail.find(query_dict, result_dict)
        gene2name = pd.DataFrame(list(result))
        gene2name.set_index('transcript_id', inplace=True)
        gene2name["type"] = "mrna"

        lnc_query_dict = dict(query_id=annot_main_id)
        lnc_result_dict = dict( _id=0, gene_name=1, gene_id=1, lncrna_id=1)
        lnc_result = lnc_annot_detail.find(lnc_query_dict, lnc_result_dict)
        if lnc_result.count() > 0:
            lnc_gene2name = pd.DataFrame(list(lnc_result))
            lnc_gene2name.rename(columns={"lncrna_id", 'transcript_id'})
            lnc_gene2name.set_index('transcript_id', inplace=True)
            lnc_gene2name["type"] = "lnc_rna"
            gene2name = gene2name.append(lnc_gene2name)

        lncnew_query_dict = dict(query_id=annot_main_id)
        lncnew_result_dict = dict(id=0, gene_name=1, gene_id=1, transcript_id=1)
        lncnew_result = lncnew_annot_detail.find(lncnew_query_dict, lncnew_result_dict)
        if lncnew_result.count() > 0:
            lncnew_gene2name = pd.DataFrame(list(lncnew_result))
            lncnew_gene2name.set_index('transcript_id', inplace=True)
            lncnew_gene2name["type"] = "lnc_rna"
            gene2name = gene2name.append(lncnew_gene2name)

    else:
        query_dict = dict(query_id=annot_main_id, is_gene=True)
        result_dict = dict( _id=0, gene_name=1, gene_id=1)
        result = annot_detail.find(query_dict, result_dict)
        gene2name = pd.DataFrame(list(result))
        gene2name.set_index('gene_id', inplace=True)
        gene2name["type"] = "mrna"

        lnc_query_dict = dict(query_id=annot_main_id, is_gene=True)
        lnc_result_dict = dict( _id=0, gene_name=1, gene_id=1)
        lnc_result = annot_detail.find(lnc_query_dict, lnc_result_dict)
        if lnc_result.count() > 0:
            lnc_gene2name = pd.DataFrame(list(lnc_result))
            lnc_gene2name.set_index('gene_id', inplace=True)
            lnc_gene2name["type"] = "lncrna"
            gene2name = gene2name.append(lnc_gene2name)

        lncnew_query_dict = dict(query_id=annot_main_id, is_gene=True)
        lncnew_result_dict = dict( _id=0, gene_name=1, gene_id=1)
        lncnew_result = annot_detail.find(lncnew_query_dict, lncnew_result_dict)
        if lncnew_result.count() > 0:
            lncnew_gene2name = pd.DataFrame(list(lncnew_result))
            lncnew_gene2name.set_index('gene_id', inplace=True)
            lncnew_gene2name["type"] = "lncrna"
            gene2name = gene2name.append(lncnew_gene2name)

    gene2name = gene2name.loc[list(target_seqs), :]
    gene2name.reset_index(inplace=True)
    output = os.path.join(dir_path, "seq_id2gene_name.txt")
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    gene2name = pd.read_table(output, header=0)
    gene2name.fillna(method="pad", axis=1, inplace=True)
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    return exp_matrix+';'+output


def export_gene_name_type(data, option_name, dir_path, bind_obj=None):
    """ 刘彬旭
    获取基因集基因name和基因集类型信息
    """
    prepare_id = data.strip()
    # prepare_main = wgcna_prepare.find_one({"_id": ObjectId(prepare_id)})
    geneset_ids = data.split(";")

    geneset = dict()
    task_id = bind_obj.sheet.option("task_id")
    geneset_type = []
    for geneset_id in geneset_ids:
        if "all" not in geneset_id.lower():
            conn = db['geneset_detail']
            conn_geneset = db['geneset']
            geneset_records1 = conn_geneset.find_one({"main_id": ObjectId(geneset_id), "task_id": task_id})
            geneset_type = geneset_records1["type"]
            geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})

            if not geneset_records:
                bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))

            if geneset_type in geneset:
                geneset[geneset_type] = geneset[geneset_type] | set(geneset_records['seq_list'])
            else:
                geneset[geneset_type] = set(geneset_records['seq_list'])

    annot_table = db['annotation_query']
    lnc_annot_table = db['known_lncrna_identify']
    lncnew_annot_table = db['new_lncrna_predict']

    try:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        lnc_annot_main = lnc_annot_table.find_one({"task_id": task_id, "status": "end"})
        lncnew_annot_main = lncnew_annot_table.find_one({"task_id": task_id, "status": "end"})
    except:
        bind_obj.set_error("cannot find annotation_query main table")
    else:
        if annot_main is None:
            annot_main = annot_table.find_one({"task_id": task_id, "type": "origin", "status": "end"})
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']


    if "main_id" not in lnc_annot_main:
        lnc_annot_main_id = lnc_annot_main['_id']
    else:
        lnc_annot_main_id = lnc_annot_main['main_id']

    if "main_id" not in lncnew_annot_main:
        lncnew_annot_main_id = lncnew_annot_main['_id']
    else:
        lncnew_annot_main_id = lncnew_annot_main['main_id']

    annot_detail = db['annotation_query_detail']
    lnc_annot_detail = db['known_lncrna_identify_detail']
    lncnew_annot_detail = db['new_lncrna_predict_detail']

    gene2name_all = pd.DataFrame()

    for gene_type, gene_set in geneset.items():
        if gene_type == "G":
            query_dict = dict(query_id=annot_main_id, is_gene=True)
            result_dict = dict( _id=0, gene_name=1, gene_id=1)
            result = annot_detail.find(query_dict, result_dict)

            gene2name = pd.DataFrame(list(result))
            gene2name.rename(columns={"gene_id": 'seq_id'}, inplace=True)
            gene2name.set_index('seq_id', inplace=True)
            gene2name["gene_type"] = "mRNA"
            gene2name_choose = gene2name[gene2name.index.isin(gene_set)]
            gene2name_all = gene2name_all.append(gene2name_choose)
            pass
        elif gene_type == "T":
            query_dict = dict(query_id=annot_main_id,)
            result_dict = dict( _id=0, gene_name=1, gene_id=1, transcript_id=1)
            result = annot_detail.find(query_dict, result_dict)
            gene2name = pd.DataFrame(list(result))
            gene2name.rename(columns={"transcript_id": 'seq_id'}, inplace=True)
            gene2name.set_index('seq_id', inplace=True)
            gene2name["gene_type"] = "mRNA"
            gene2name_choose = gene2name[gene2name.index.isin(gene_set)]
            gene2name_all = gene2name_all.append(gene2name_choose)

        elif gene_type == "LG":
            lnc_query_dict = dict(query_id=annot_main_id, is_gene=True)
            lnc_result_dict = dict( _id=0, gene_name=1, gene_id=1)
            lnc_result = annot_detail.find(lnc_query_dict, lnc_result_dict)
            lnc_gene2name = pd.DataFrame(list(lnc_result))
            lnc_gene2name.rename(columns={"gene_id": 'seq_id'}, inplace=True)
            # lnc_gene2name.to_csv("~/aaa" + "test.tmp", sep='\t', header=True, index=False)
            lnc_gene2name.set_index('seq_id', inplace=True)
            lnc_gene2name["gene_type"] = "lncRNA"

            lncnew_query_dict = dict(query_id=annot_main_id, is_gene=True)
            lncnew_result_dict = dict( _id=0, gene_name=1, gene_id=1)
            lncnew_result = annot_detail.find(lncnew_query_dict, lncnew_result_dict)
            lncnew_gene2name = pd.DataFrame(list(lncnew_result))
            lncnew_gene2name.rename(columns={"gene_id": 'seq_id'}, inplace=True)
            lncnew_gene2name.set_index('seq_id', inplace=True)
            lncnew_gene2name["gene_type"] = "lncRNA"

            lnc_gene2name_choose = lnc_gene2name[lnc_gene2name.index.isin(gene_set)]
            lncnew_gene2name_choose = lncnew_gene2name[lncnew_gene2name.index.isin(gene_set)]

            gene2name_all = gene2name_all.append(lnc_gene2name_choose)
            gene2name_all = gene2name_all.append(lncnew_gene2name_choose)
        elif gene_type == "LT":
            lnc_query_dict = dict(query_id=annot_main_id)
            lnc_result_dict = dict( _id=0, gene_name=1, gene_id=1, lncrna_id=1)
            lnc_result = lnc_annot_detail.find(lnc_query_dict, lnc_result_dict)
            lnc_gene2name = pd.DataFrame(list(lnc_result))
            lnc_gene2name.rename(columns={"lncrna_id": 'seq_id'}, inplace=True)
            lnc_gene2name.set_index('seq_id', inplace=True)
            lnc_gene2name["gene_type"] = "lncRNA"
            lncnew_query_dict = dict(query_id=annot_main_id)
            lncnew_result = lncnew_annot_detail.find(lncnew_query_dict, lncnew_result_dict)
            lncnew_gene2name = pd.DataFrame(list(lncnew_result))
            lncnew_gene2name.rename(columns={"transcript_id": 'seq_id'}, inplace=True)
            lncnew_gene2name.set_index('seq_id', inplace=True)
            lncnew_gene2name["gene_type"] = "lncRNA"

            lnc_gene2name_choose = lnc_gene2name[lnc_gene2name.index.isin(gene_set)]
            lncnew_gene2name_choose = lncnew_gene2name[lncnew_gene2name.index.isin(gene_set)]

            gene2name_all = gene2name_all.append(lnc_gene2name_choose)
            gene2name_all = gene2name_all.append(lncnew_gene2name_choose)

    gene2name_all.reset_index(inplace=True)
    output = os.path.join(dir_path, "seq_id2gene_name.txt")
    gene2name_all.to_csv(output, sep='\t', header=True, index=False)
    return output

def export_gene_name_type2(data, option_name, dir_path, bind_obj=None):
    """ 刘彬旭
    获取基因集基因name和基因集类型信息
    """
    prepare_id = data.strip()
    # prepare_main = wgcna_prepare.find_one({"_id": ObjectId(prepare_id)})
    geneset_ids = data.split(";")
    mrnaset_id = geneset_ids[0]
    lncrnaset_id = geneset_ids[1]
    gene_type = bind_obj.sheet.option("exp_level")

    geneset = dict()
    task_id = bind_obj.sheet.option("task_id")
    geneset_type = []
    conn = db['geneset_detail']
    conn_geneset = db['geneset']
    if mrnaset_id != "" and mrnaset_id.lower() != "all" and mrnaset_id.lower() != "refall":
        mrnaset_records = conn.find_one({"geneset_id": ObjectId(mrnaset_id)})
        mrnaset_list = mrnaset_records['seq_list']
        mrnaset_records1 = conn_geneset.find_one({"main_id": ObjectId(mrnaset_id), "task_id": task_id})
        mrnaset_type = mrnaset_records1["type"]
    else:
        mrnaset_list = []
    if lncrnaset_id != "" and lncrnaset_id.lower() != "all" and mrnaset_id.lower() != "refall":
        lncrnaset_records = conn.find_one({"geneset_id": ObjectId(lncrnaset_id)})
        lncrnaset_list = lncrnaset_records['seq_list']
        lncrnaset_records1 = conn_geneset.find_one({"main_id": ObjectId(lncrnaset_id), "task_id": task_id})
        lncrnaset_type = lncrnaset_records1["type"]
    else:
        lncrnaset_list = []

    genes_col = db['genes']
    genes_detail_col = db['genes_detail']
    genes_main = genes_col.find_one({"task_id": task_id})

    output = os.path.join(dir_path, "seq_id2gene_name.txt")
    with open(output, 'w') as id_name_f:
        if gene_type == "G":
            id_name_f.write("seq_id\tgene_type\tgene_name\n")
            all_mrna_detail = genes_detail_col.find({"genes_id": genes_main["main_id"], "gene_type" : "mRNA", "level" : "G"})
            if mrnaset_id.lower() in ["all", "refall"]:
                for a_gene in all_mrna_detail:
                    if a_gene["gene_name"][0] in ["-", ""]:
                        gene_name = a_gene["gene_id"]
                    else:
                        gene_name = a_gene["gene_name"][0]
                    id_name_f.write("{}\t{}\t{}\n".format(
                        a_gene["gene_id"],
                        a_gene["gene_type"],
                        gene_name
                    ))
            elif mrnaset_id.lower() != "":
                for a_gene in all_mrna_detail:
                    if a_gene["gene_id"] in mrnaset_list:
                        if a_gene["gene_name"][0] in ["-", ""]:
                            gene_name = a_gene["gene_id"]
                        else:
                            gene_name = a_gene["gene_name"][0]
                        id_name_f.write("{}\t{}\t{}\n".format(
                            a_gene["gene_id"],
                            a_gene["gene_type"],
                            gene_name
                        ))
            all_lncrna_detail = genes_detail_col.find({"genes_id": genes_main["main_id"], "gene_type" : "lncRNA", "level" : "G"})
            if lncrnaset_id.lower() in ["all", "refall"]:
                for a_gene in all_lncrna_detail:
                    if a_gene["gene_name"][0] in ["-", ""]:
                        gene_name = a_gene["gene_id"]
                    else:
                        gene_name = a_gene["gene_name"][0]
                    id_name_f.write("{}\t{}\t{}\n".format(
                        a_gene["gene_id"],
                        a_gene["gene_type"],
                        gene_name
                    ))
            elif lncrnaset_id.lower() != "":
                for a_gene in all_lncrna_detail:
                    if a_gene["gene_id"] in lncrnaset_list:
                        if a_gene["gene_name"][0] in ["-", ""]:
                            gene_name = a_gene["gene_id"]
                        else:
                            gene_name = a_gene["gene_name"][0]
                        id_name_f.write("{}\t{}\t{}\n".format(
                            a_gene["gene_id"],
                            a_gene["gene_type"],
                            gene_name
                        ))
        elif gene_type == "T":
            id_name_f.write("seq_id\tgene_id\tgene_type\tgene_name\n")
            all_mrna_detail = genes_detail_col.find({"genes_id": genes_main["main_id"], "gene_type" : "mRNA", "level" : "T"})
            if mrnaset_id.lower() in ["all", "refall"]:
                for a_gene in all_mrna_detail:
                    if a_gene["gene_name"][0] in ["-", ""]:
                        gene_name = a_gene["gene_id"]
                    else:
                        gene_name = a_gene["gene_name"][0]
                    id_name_f.write("{}\t{}\t{}\t{}\n".format(
                        a_gene["trans_id"],
                        a_gene["gene_id"],
                        a_gene["gene_type"],
                        gene_name
                    ))
            elif mrnaset_id.lower() != "":
                for a_gene in all_mrna_detail:
                    if a_gene["trans_id"] in mrnaset_list:
                        if a_gene["gene_name"][0] in ["-", ""]:
                            gene_name = a_gene["gene_id"]
                        else:
                            gene_name = a_gene["gene_name"][0]
                        id_name_f.write("{}\t{}\t{}\t{}\n".format(
                            a_gene["trans_id"],
                            a_gene["gene_id"],
                            a_gene["gene_type"],
                            gene_name
                        ))
            all_lncrna_detail = genes_detail_col.find({"genes_id": genes_main["main_id"], "gene_type" : "lncRNA", "level" : "T"})
            if lncrnaset_id.lower() in ["all", "refall"]:
                for a_gene in all_lncrna_detail:
                    if a_gene["gene_name"][0] in ["-", ""]:
                        gene_name = a_gene["gene_id"]
                    else:
                        gene_name = a_gene["gene_name"][0]
                    id_name_f.write("{}\t{}\t{}\t{}\n".format(
                        a_gene["trans_id"],
                        a_gene["gene_id"],
                        a_gene["gene_type"],
                        a_gene["gene_name"][0]
                    ))
            elif lncrnaset_id.lower() != "":
                for a_gene in all_lncrna_detail:
                    if a_gene["trans_id"] in lncrnaset_list:
                        if a_gene["gene_name"][0] in ["-", ""]:
                            gene_name = a_gene["gene_id"]
                        else:
                            gene_name = a_gene["gene_name"][0]
                        id_name_f.write("{}\t{}\t{}\t{}\n".format(
                            a_gene["trans_id"],
                            a_gene["gene_id"],
                            a_gene["gene_type"],
                            gene_name
                        ))

    return output

def export_genename_des(data, option_name, dir_path, bind_obj=None):
    """ 刘彬旭
    获取基因集基因name和描述信息
    """
    task_id = bind_obj.sheet.option("task_id")
    gene_level = bind_obj.sheet.option("type")
    bind_obj.logger.info("**** data")
    print "data is {}".format(data)
    geneset_ids = data.split(";")
    bind_obj.logger.info("{}".format(geneset_ids))
    conn = db['geneset_detail']
    genes_col = db['genes']
    geneset = set()
    for geneset_id in geneset_ids:
        bind_obj.logger.info("{}".format(geneset_id))
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        geneset_list = geneset_records['seq_list']
        bind_obj.logger.info("{}".format(geneset_id))
        geneset = geneset | set(geneset_list)

    # bind_obj.logger.info("gene set is {}".format(geneset))

    genes_detail_col = db['genes_detail']

    genes_main = genes_col.find_one({"task_id": task_id})
    output = os.path.join(dir_path, "geneset_annot.xls")
    with open(output, 'w') as geneannot_f:
        geneannot_f.write("gene_id\ttranscript_id\tis_gene\tgene_name\tlength\tdescription\n")
        for gene in geneset:
            if gene_level == "T":
                genes_detail_dict = genes_detail_col.find({"genes_id": genes_main["main_id"], "trans_id" : gene, "level" : "T"})
            else:
                bind_obj.logger.info("gene set is {} {}".format(genes_main["main_id"], gene))
                print genes_main["main_id"], gene
                genes_detail_dict = genes_detail_col.find({"genes_id": genes_main["main_id"], "gene_id" : gene, "level" : "T"})
            for a_gene in  genes_detail_dict:
                # bind_obj.logger.info("gene set is {}".format(a_gene))
                # bind_obj.logger.info("gene set is {}".format(a_gene["description"]))
                geneannot_f.write("\t".join([a_gene["gene_id"],
                                             a_gene["trans_id"],
                                             "yes",
                                             a_gene["gene_name"][0],
                                             str(a_gene["gene_length"]),
                                             str(a_gene["description"][0])
                ]) + "\n")

    return output


# added  by gdq for wgcna
def export_wgcna_relate_input(data, option_name, dir_path, bind_obj=None):
    module_id = data.strip()
    # export eigengenes
    eigengenes = db['wgcna_module_eigengenes_detail']
    eigengenes_found = eigengenes.find({"module_id": ObjectId(module_id)}, {"_id": 0, "module_id":0})
    eigengenes_pd = pd.DataFrame(list(eigengenes_found))
    eigengenes_pd.set_index("module", inplace=True)
    eigengenes_path = os.path.join(dir_path, "module_eigengenes.xls")
    eigengenes_pd.to_csv(eigengenes_path, sep='\t', header=True, index=True)
    # export exp matrix
    prepare_id = db['wgcna_module'].find_one({"main_id": ObjectId(module_id)})['wgcna_prepare_id']
    exp_matrix = db['wgcna_prepare'].find_one({"main_id": ObjectId(prepare_id)})['exp_matrix']
    # export each gene's module and gene_id/gene_name info
    exp_level = bind_obj.sheet.option('exp_level')
    membership = db['wgcna_module_membership_detail']
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
    conn = db['specimen']
    result=conn.find({"task_id": data, "about_qc": "after"})
    sorted_result = sorted(result, key=lambda k: k['_id'])
    output = os.path.join(dir_path, option_name)
    with open(output,'w') as f:
        for i in sorted_result:
            f.write(i['bam_path'] + "\n")
    return output

#---------------------从蛋白复制过来的------------------------

def export_compare_exp_fc(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析accession_id 和 fc
    '''

    diff_id = bind_obj.sheet.option('diff_id')
    compare_group = bind_obj.sheet.option('compare_group')
    target_cols = OrderedDict(seq_id=1, log2fc=1, _id=0)

    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id),"compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    output = os.path.join(dir_path, option_name)
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export expression matrix')
    return output

def export_go_enrich_matrix(data, option_name, dir_path, bind_obj=None):

    go_enrich_id = bind_obj.sheet.option('enrich_id')
    go_type = bind_obj.sheet.option('go_type')

    target_cols = OrderedDict(go_id=1, go_type=1, discription=1, p_corrected=1, p_uncorrected=1, seq_list=1, depth=1, _id=0)
    bind_obj.logger.debug("导出GO富集表")
    # get geneset
    conn = db['geneset_go_enrich_detail']
    if go_type == 'ALL' or go_type == 'All':
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id)}, target_cols)
    else:
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id), "go_type":go_type}, target_cols)
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
    conn = db['geneset_kegg_enrich_detail']
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
    collection = db['geneset_detail']
    main_collection = db['geneset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:{}在geneset中未找到！".format(ObjectId(data)))
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
        print 'file can not find'
    return to_path

# added for transcription factor analysis  --------------------------------------------
def get_all_pep_seq(data, option_name, dir_path, bind_obj=None):
    pep_db_path = data.strip()
    result_path = os.path.join(dir_path, "all_pep.fa")
    if re.match(r'^\w+://\S+/.+$', pep_db_path) or re.match(r'/mnt/ilustre', pep_db_path):
        transfer = MultiFileTransfer()
        transfer.add_download(pep_db_path, bind_obj.work_dir + "/")
        transfer.perform()
        #download_s3_file(pep_db_path, os.path.join(bind_obj.work_dir, "pep.db"))
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
    annot_table = db['annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", 'status': 'end'})
    if not annot_main:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    if not annot_main:
        task_table = db['task']
        task_main = task_table.find_one({"task_id": task_id})
        if os.path.exists(task_main['annot']):
            gene_annot = os.path.join(dir_path, "seq_annot.xls")
            os.system("cut -f 1-4 " + task_main['annot'] +" |awk '{print $2\"\\t\"$1\"\\t\"$4}' > " + gene_annot)
            return gene_annot
            # os.system("cut -f 1-4 {} |awk '{print $2\"\t\"$1\"\t\"$4}' > {} ".format(task_main['annot'], gene_annot))
        else:
            bind_obj.set_error("Not Found in annotation_query by query {}".format(task_id))
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['annotation_query_detail']
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
    os.system(r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot


# added by gdq for tfbs predict
def export_geneid2tfid_file(data, option_name, dir_path, bind_obj=None):
    conn = db['tf_predict_detail']
    result = conn.find({"tf_predict_id": ObjectId(data.strip())}, {"_id": 0, 'gene_id': 1, "blast_hit": 1})
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("gene_id", inplace=True)
    result_file = os.path.join(dir_path, "geneid2tfid.txt")
    result_pd.to_csv(result_file, sep='\t', header=True, index=True)
    return result_file


def export_predict_result(data, option_name, dir_path, bind_obj=None):
    conn = db['tfbs_predict_detail']
    target_cols = {
        "_id":0, "tf_geneid":1, "gene_name_tf":1, "target_id":1, "gene_name_target":1,
        "corr":1, "corr_pvalue":1, "corr_padjust":1, "p-value":1, "q-value":1
    }
    result = conn.find({"tfbs_predict_id": ObjectId(data.strip()),}, target_cols,)
    result_pd = pd.DataFrame(list(result))
    result_file = os.path.join(dir_path, "tfbs_predict.xls")
    result_pd.to_csv(result_file, sep='\t', header=True, index=False)
    return result_file

# ----------------------------------------------------------------------------------------------------------------------

def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} = {}'.format(k, v))

def export_bam_list(data, option_name, dir_path, bind_obj=None):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    results = db['specimen'].find({'task_id': data, 'about_qc': 'after'})
    output = os.path.join(dir_path, 'bam.list')
    open(output, 'w').writelines(sorted(['{}\n'.format(i['bam_path']) for i in results]))
    return output

def export_rmats_group_table(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    if isinstance(data, types.StringTypes):
        group_dict = json.loads(data)
    elif isinstance(data, types.DictType):
        group_dict = data
    lines = ['#sample\tgroup\n']
    for group, samples in group_dict.items():
        lines.extend(['{}\t{}\n'.format(sample, group) for sample in sorted(samples)])
    group_table = os.path.join(dir_path, 'group.txt')
    open(group_table, 'w').writelines(lines)
    return group_table

def export_rmats_control_table(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    control_table = os.path.join(dir_path, 'control.txt')
    open(control_table, 'w').writelines(['#control\tother\n', '{}\t{}\n'.format(*reversed(data.split('|')))])
    return control_table

def export_rmats_detail_loc2name(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    lines = list()
    if 'task_id' in bind_obj.sheet.data:
        task_id = bind_obj.sheet.data['task_id']
    else:
        task_id = "_".join(bind_obj.sheet.data['id'].split("_")[0:-2])
    for splicing_id in data.split(','):
        rmats_info = db['splicing_rmats'].find_one(
            {'main_id': ObjectId(splicing_id), 'task_id': task_id}
        )
        lines.append('{}\t{}\n'.format(
            os.path.join(rmats_info['result_dir'], 'all_events_detail_big_table.txt'), rmats_info['compare_plan']
        ))
    output = os.path.join(dir_path, 'rmats_detail.list')
    open(output, 'w').writelines(lines)
    return output

def export_lot_of_lncrna_list(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    geneset_detail_info = db['geneset_detail'].find_one({'geneset_id': ObjectId(data)})
    output = os.path.join(dir_path, 'keep.m.list')
    open(output, 'w').writelines('{}\n'.format(i) for i in geneset_detail_info['seq_list'])
    return output

def export_single_lncrna_list(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    output = os.path.join(dir_path, 'keep.s.list')
    open(output, 'w').write('{}\n'.format(data))
    return output

def export_target_fa(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    target = os.path.join(bind_obj.config.SOFTWARE_DIR, 'database/Genome_DB_finish', data, 'lncrna.fa')
    output = os.path.join(dir_path, 'target.fa')
    bind_obj.logger.debug('start copying {} to {}'.format(target, output))
    shutil.copy(target, output)
    return output

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
    with open(dir_path+'/group_info.txt', 'w') as f:
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
            conn = db['geneset_detail']
            geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
            if not geneset_records:
                bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
            geneset.update(geneset_records['seq_list'])

    # get all exp matrix
    conn = db['exp_detail']
    exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: {}'.format(exp_id))
    exp_matrix = exp_matrix.set_index('seq_id')
    if geneset:
        exp_matrix = exp_matrix.loc[geneset, :]
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
    bind_obj.logger.debug('=========== %s ++++++++++++++' % ','.join(data))
    geneset_ids = [item.split(',') for item in data.split(";")]
    geneset = {}
    for gene_type, geneset_id in geneset_ids:
        if "all" not in geneset_id.lower():
            conn = db['geneset_detail']
            geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
            if not geneset_records:
                bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
            for i in geneset_records['seq_list']:
                geneset[i] = 'L' if gene_type == 'l_geneset_id' else 'G'

    output = os.path.join(dir_path, 'genes_info.json')
    with open(output, 'w') as in_handler:
        json.dump(geneset, in_handler, indent=4)

    return output
