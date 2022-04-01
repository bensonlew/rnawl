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
    # data = json.loads(data)
    # geneset_id_list = data['geneset_id'].split(",")
    geneset_id_list = data.split(",")
    source = bind_obj.sheet.option('source')
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


def export_do_anno_genesets(data, option_name, dir_path, bind_obj=None):
    # {'c1': data.c1, 'c2': data.c2, 'c3': data.c3, 'task_id': data.task_id,
    #                                   'type': 'origin'}
    bind_obj.logger.info("data is {}".format(data))
    task_id = bind_obj.sheet.option('task_id')
    genesets = data.split(",")
    source = bind_obj.sheet.option('source')

    if False:
        geneset_id = genesets[0]
        genesets_coll = db['sg_geneset']
        my_result = genesets_coll.find_one({'main_id': ObjectId(geneset_id)})
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
        do_files = list()
        for n, gi in enumerate(genesets):
            genesets_coll = db['sg_geneset']
            genesets_detail_coll = db['sg_geneset_detail']
            my_result = genesets_coll.find_one({'main_id': ObjectId(gi)})

            if not my_result:
                bind_obj.set_error("意外错误，geneset_id:%s在sg_geneset中未找到！", variables=(ObjectId(gi)))

            results = genesets_detail_coll.find_one({"geneset_id": ObjectId(gi)})
            genes = results["seq_list"]

            do_file = dir_path + '/' +  my_result["name"] + '.do.list'
            do_files.append(do_file)
            exp_coll = db['sg_exp']
            exp_result = exp_coll.find_one({'task_id': bind_obj.sheet.option('task_id'), 'is_rmbe':False, 'level': "G"})
            exp_id = exp_result["main_id"]
            exp_detail_coll = db['sg_exp_detail']
            target_cols = OrderedDict(gene_id=1, do=1, _id=0)
            exp_records = exp_detail_coll.find({"exp_id": ObjectId(exp_id), "gene_id": {"$in": genes}}, target_cols)
            with open(do_file, 'w') as fo:
                for rec in exp_records:
                    if "do" in rec:
                        if rec["do"] == "":
                            pass
                        else:
                            do_list = rec["do"].split("; ")
                            do_clean = [x.split("(")[0] for x in do_list]
                            fo.write(rec['gene_id'] + '\t' + ";".join(do_clean) + '\n')
        return ",".join(do_files)

def export_do_list(data, option_name, dir_path, bind_obj=None):
    '''
    get all gene do annot list
    '''
    bind_obj.logger.info("data is {}".format(data))
    task_id = bind_obj.sheet.option('task_id')
    exp_coll = db['sg_exp']
    exp_result = exp_coll.find_one({'task_id': bind_obj.sheet.option('task_id'), 'level': "G", 'is_rmbe': False})
    exp_id = exp_result["main_id"]
    exp_detail_coll = db['sg_exp_detail']
    target_cols = OrderedDict(gene_id=1, do=1, _id=0)
    exp_records = exp_detail_coll.find({"exp_id": ObjectId(exp_id), }, target_cols)
    do_files = dir_path + '/all_do.list'
    with open(do_files, 'w') as fo:
        for rec in exp_records:
            if rec["do"] == "":
                pass
            else:
                do_list = rec["do"].split("; ")
                do_clean = [x.split("(")[0] for x in do_list]
                fo.write(rec['gene_id'] + '\t' + ";".join(do_clean) + '\n')

    return do_files


def export_reactome_annot(data, option_name, dir_path, bind_obj=None):
    '''
    get all gene do annot list
    '''
    bind_obj.logger.info("data is {}".format(data))
    task_id = bind_obj.sheet.option('task_id')
    exp_coll = db['sg_exp']
    exp_result = exp_coll.find_one({'task_id': bind_obj.sheet.option('task_id'), 'level': "G", 'is_rmbe': False})
    exp_id = exp_result["main_id"]
    exp_detail_coll = db['sg_exp_detail']
    target_cols = OrderedDict(gene_id=1, reactome_link=1, _id=0)
    exp_records = exp_detail_coll.find({"exp_id": ObjectId(exp_id), }, target_cols)
    reac_files = dir_path + '/all_reactome.list'
    with open(reac_files, 'w') as fo:
        for rec in exp_records:
            if rec["reactome_link"] == "":
                pass
            else:
                fo.write(rec['gene_id'] + '\t' + rec["reactome_link"] + '\n')

    return reac_files



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
