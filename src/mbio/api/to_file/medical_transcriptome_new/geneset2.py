# -*- coding: utf-8 -*-
# __author__ = 'sanger'

from __future__ import division
import os
from itertools import chain

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
# from biocluster.api.file.lib.s3 import S3TransferManager
# from boto.s3.bucket import Bucket
import shutil
from biocluster.api.file.lib.transfer import MultiFileTransfer
import sys

project_type = 'medical_transcriptome'
project_type1 = 'ref_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
# msi_db = Config().get_mongo_client(mtype=project_type, ref=True)[Config().get_mongo_dbname(project_type,ref=True)]
ref_rna_v2_db = Config().get_mongo_client(mtype=project_type1)[Config().get_mongo_dbname(project_type1)]

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
                bind_obj.set_error("意外错误，geneset_id:%s在sg_geneset中未找到！".format(ObjectId(geneset_id)))
                # bind_obj.set_error("意外错误，geneset_id:%s在sg_geneset中未找到！", variables=(ObjectId(geneset_id)), code="55700108")
            geneset_name = my_result["name"]
            results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
            if not results:
                bind_obj.set_error("意外错误，基因集详情在sg_geneset_detail中未找到！请核实基因集来源并联系技术")
            seq_list = results["seq_list"]
            if 'regulate_list' not in results:
                bind_obj.set_error('geneset_detail table must contain "regulate_list" field when provide one geneset id')
                # bind_obj.set_error(
                #     'geneset_detail table must contain "regulate_list" field when provide one geneset id')
            regulate_list = results["regulate_list"]
            up_list, down_list = [], []
            for seq_id, regulate in zip(seq_list, regulate_list):
                if regulate == 'up':
                    up_list.append(seq_id)
                else:
                    down_list.append(seq_id)
            out_handler.write(geneset_name + '_up\t' + ','.join(up_list) + '\n')
            out_handler.write(geneset_name + '_down\t' + ','.join(down_list) + '\n')
        else:
            for n, gi in enumerate(geneset_id_list):
                my_result = main_collection.find_one({'main_id': ObjectId(gi)})
                if not my_result:
                    bind_obj.set_error("意外错误，geneset_id:%s在sg_geneset中未找到！", variables=(ObjectId(gi)), code="55700110")
                out_handler.write(my_result["name"] + "\t")
                results = collection.find_one({"geneset_id": ObjectId(gi)})
                out_handler.write(",".join(results["seq_list"]) + "\n")

    return multi_geneset_path


# 通过基因集id到sg_geneset去获取task_id，type对应的信息，然后到sg_annotation_kegg去查找注释信息，然后导表；页面上的选择限制2个geneset_id的类型
# 肯定是一样，所以任意选择一个geneset_id过来就可以，所以接口那里随便选择了一个
# 这里的one_record是多余的，直接find的结果就可以判断，既然写了就写了吧






# 根据"add_info":geneset_info['task_id'] + "\t" + data.geneset_type，就是sg_geneset主表的task_id和页面传过来的是annotation还是annotation1
# 也就是是origin还是latest









def export_genes_detail_medical(data, option_name, dir_path, bind_obj=None):
    task_id, exp_level = data['task_id'], data['exp_level']
    col = db["sg_exp"]
    main_id = col.find_one({"task_id": task_id, "exp_level": exp_level})['main_id']

    col = db["sg_exp_detail"]
    if exp_level == "G":
        choose_fields = {'gene_id': 1, 'gene_name': 1, 'description': 1, '_id': 0}
    else:
        choose_fields = {'transcript_id': 1, 'gene_name': 1, 'description': 1, '_id': 0}

    results = col.find({"exp_id": ObjectId(main_id)}, choose_fields)
    # Probe Set ID    Gene Symbol     Gene Title
    out_list = []
    line_demo = '{gene_id}\t{gene_name}\t{description}\n'
    seq_id_set = set()
    out_file = os.path.join(dir_path, 'gsea.chip')
    try:
        geneset_source = bind_obj.sheet.option('geneset_source')
    except:
        geneset_source = ""


    with open(out_file, 'w') as out_handler:
        out_handler.write('gene_id\tgene_name\tdescription\n')
        for dic in results:
            if "gene_name" in dic:
                pass
            else:
                continue

            if geneset_source == 'msigdb' and dic["gene_name"] in [None, "None", "", "-", 'nan']:
                continue
            if "description" in dic:
                if dic["description"] in [None, "None", "", 'nan']:
                    dic["description"] = "-"
            else:
                dic["description"] = "-"
            if exp_level == "T":
                gene_id = dic['transcript_id']
            else:
                gene_id = dic['gene_id']
            if gene_id in seq_id_set:
                continue
            seq_id_set.add(gene_id)
            # bind_obj.logger.info("data is {}".format(dic))
            if exp_level == "G":
                try:
                    out_handler.write(line_demo.format(gene_id=dic["gene_id"], gene_name=dic["gene_name"].upper(), description=dic["description"]))
                except:
                    print dic["gene_name"]
            if exp_level == "T":
                out_handler.write(line_demo.format(gene_id=dic["transcript_id"], gene_name=dic["gene_name"].upper(), description=dic["description"]))

    return out_file







def export_do_anno_genesets(data, option_name, dir_path, bind_obj=None):
    # {'c1': data.c1, 'c2': data.c2, 'c3': data.c3, 'task_id': data.task_id,
    #                                   'type': 'origin'}
    bind_obj.logger.info("data is {}".format(data))
    task_id = bind_obj.sheet.option('task_id')
    # anno_level = bind_obj.sheet.option('level')
    genesets = data.split(",")




    go_list = {i.strip() for i in genesets.strip().split(',')}
        '''
        out_file = os.path.join(dir_path, 'go.list')
        with open(out_file, 'w') as out_handler:
            out_handler.write(
                ''.join(['{}\t{}\n'.format(dic['gene_id'], dic['gos_list']) for dic in results])
            )
        return out_file
    else:
        main_id = db['sg_annotation_kegg'].find_one({'task_id': task_id, 'type': "origin"})['main_id']
        query_dic = {'kegg_id': main_id, 'anno_type': anno_level}
        if subtype == "defined":
            pass
        else:
            if subtype2 not in ['all', 'All', ""]:
                if len(subtype2.split(",")) > 1:
                    subtype2_par = subtype2.split(",")
                    subtype2_list = [s.replace("_", ",") for s in subtype2_par]
                    query_dic['second_category'] = {"$in": subtype2_list}
                else:
                    query_dic['second_category'] = subtype2.replace("_", ",")

        query_dic['seq_type'] = "all"


        results = db['sg_annotation_kegg_level'].find(query_dic, {'seq_list': 1, 'pathway_id': 1, 'pathway_definition': 1, '_id': 0})
        kegg_list = None
        if genesets not in ["all", ""]:
            kegg_list = {i.strip() for i in genesets.strip().split(',')}
        out_file = os.path.join(dir_path, 'kegg.gmt')
        with open(out_file, 'w') as out_handler:
            if kegg_list:
                out_handler.write(
                    ''.join('{}\t{}\t{}\n'.format(
                        dic['pathway_id'], dic['pathway_definition'], dic['seq_list'].replace(';', '\t')) for dic in results if dic['pathway_id'] in kegg_list)
                )
            else:
                out_handler.write(
                    ''.join(['{}\t{}\t{}\n'.format(
                        dic['pathway_id'], dic['pathway_definition'], dic['seq_list'].replace(';', '\t')) for dic in results])
                )
        return out_file


def export_msigdb_genesets(data, option_name, dir_path, bind_obj=None):
    # {'c1': data.c1, 'c2': data.c2, 'c3': data.c3}
    bind_obj.logger.info("data is {}".format(data))
    extr_dic = {}
    species= data.split(";")[0]
    extr_dic['organism'] = species
    c1 = data.split(";")[1]
    extr_dic['c1'] = c1
    c2 = data.split(";")[2]
    if c2.count('|') == 1:
        c2, tmp_c3 = c2.split('|')
        extr_dic['c2'] = c2
        extr_dic['c3'] = tmp_c3
    else:
        extr_dic['c2'] = c2

    results = ref_rna_v2_db['sg_msigdb'].find(extr_dic, {'_id': 0})
    ref_link = re.compile('<[^>]*>')
    res_dict = {dic['name']: '{}\t{}\t{}\n'.format(
        dic["name"],
        re.sub(ref_link, '', dic['brief_description']),
        dic["gene_symbols"].replace(" ", "\t")) for dic in results if 'name' in dic}
    c3 = data.split(";")[3]
    out_file = os.path.join(dir_path, 'gsea.gmt')
    with open(out_file, 'w') as out_handler:
        if c3 == 'all' or c3 =='':
            out_handler.write(''.join(res_dict.values()))
        else:
            keys = {i.strip() for i in c3.split(',')}
            out_handler.write(''.join(v for k, v in res_dict.items() if k in keys))
    return out_file

def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} - {}'.format(k, v))
    else:
        pass
