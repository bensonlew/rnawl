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
project_type1 = 'ref_rna'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
# msi_db = Config().get_mongo_client(mtype=project_type, ref=True)[Config().get_mongo_dbname(project_type,ref=True)]
medical_transcriptome_db = Config().get_mongo_client(mtype=project_type1, db_version=0, ref=True)[Config().get_mongo_dbname(project_type1, db_version=0, ref=True)]

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
def export_kegg_table(data, option_name, dir_path, bind_obj=None):
    kegg_path = os.path.join(dir_path, 'gene_kegg_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg_path))
    geneset_collection = db["sg_geneset"]
    bind_obj.logger.debug(data)
    geneset_result = geneset_collection.find_one({"main_id": ObjectId(data)})
    task_id = geneset_result["task_id"]
    bind_obj.logger.debug("ttttttt")
    bind_obj.logger.debug(task_id)
    geneset_type = geneset_result["type"]
    my_result = db["sg_annotation_kegg"].find(
        {"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\n')
        for main_table in my_result:
            kegg_id = main_table["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not my_result:
                bind_obj.set_error("意外错误，annotation_kegg_id:%s在sg_annotation_kegg中未找到！", variables=(kegg_id), code="55700111")
            results = db['sg_annotation_kegg_table'].find({'kegg_id': kegg_id, 'anno_type': geneset_type})
            one_record = db['sg_annotation_kegg_table'].find_one({'kegg_id': kegg_id, 'anno_type': geneset_type})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:%s在sg_annotation_kegg_table中未找到！", variables=(ObjectId(kegg_id)), code="55700112")
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['ko_id'] + result['transcript_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'],
                                                      result['hyperlink'], result['paths']))
    return kegg_path


def export_kegg_level_table(data, option_name, dir_path, bind_obj=None):
    kegg__level_path = os.path.join(dir_path, 'gene_kegg_level_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg__level_path))
    geneset_collection = db["sg_geneset"]
    bind_obj.logger.debug(data)
    geneset_result = geneset_collection.find_one({"main_id": ObjectId(data)})
    task_id = geneset_result["task_id"]
    bind_obj.logger.debug(task_id)
    my_result = db["sg_annotation_kegg"].find(
        {"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    with open(kegg__level_path, 'wb') as w:
        w.write(
            'Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
        for i in my_result:
            kegg_id = i["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not kegg_id:
                bind_obj.set_error("意外错误，annotation_kegg_id:%s在sg_annotation_kegg中未找到！", variables=(kegg_id), code="55700113")
            results = db["sg_annotation_kegg_level"].find(
                {"kegg_id": kegg_id, "seq_type": "all", "anno_type": bind_obj.sheet.option("geneset_type")})
            one_record = db['sg_annotation_kegg_level'].find_one(
                {'kegg_id': kegg_id, "seq_type": "all", 'anno_type': bind_obj.sheet.option("geneset_type")})
            if not one_record:
                results = db["sg_annotation_kegg_level"].find(
                    {"kegg_id": kegg_id, "seq_type": "ref", "anno_type": bind_obj.sheet.option("geneset_type")})
                one_record = db['sg_annotation_kegg_level'].find_one(
                    {'kegg_id': kegg_id, "seq_type": "ref", 'anno_type': bind_obj.sheet.option("geneset_type")})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:%s在sg_annotation_kegg_level中未找到！", variables=(ObjectId(kegg_id)), code="55700114")
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['pathway_id'] + result['graph_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '',
                                                                          result['number_of_seqs'],
                                                                          result['pathway_definition'],
                                                                          result['first_category'], result['anno_type'],
                                                                          result['hyperlink'], result['seq_list'], '',
                                                                          result['second_category']))
    return kegg__level_path


# 根据"add_info":geneset_info['task_id'] + "\t" + data.geneset_type，就是sg_geneset主表的task_id和页面传过来的是annotation还是annotation1
# 也就是是origin还是latest
def export_add_info(data, option_name, dir_path, bind_obj=None):
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出add_info信息")
    col = db["sg_annotation_kegg"]
    result = col.find_one({"task_id": task_id, "type": bind_obj.sheet.option('type')})
    insert_id = result["main_id"]
    col = db["sg_annotation_kegg_level"]
    results = col.find({"kegg_id": insert_id, "anno_type": anno_type})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info


def export_geneset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)

    exp_id, exp_level, is_rmbe,geneset_id = data.split(';')
    conn = db['sg_exp_detail']
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    samples_raw = list()
    if is_rmbe == 'false':
        for each in group_dict:
            samples += group_dict[each]
    if is_rmbe == 'true':
        for each in group_dict:
            samples_raw += group_dict[each]
            samples += ['{}_batch'.format(x) for x in group_dict[each]]
    if exp_level == 'G':
        target_cols = OrderedDict(gene_id=1, _id=0)
    if exp_level == 'T':
        target_cols = OrderedDict(transcript_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    if 'kind' not in bind_obj.sheet.options() or bind_obj.sheet.options('kind') =='all':
        exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    else:
        if bind_obj.sheet.option('kind') in ["ref", 'new']:
            exp_records = conn.find({"exp_id": ObjectId(exp_id), 'kind': bind_obj.sheet.option('kind')}, target_cols)
    df = pd.DataFrame(list(exp_records))
    if exp_level == 'T':
        df.rename(columns={'transcript_id': 'seq_id'}, inplace=True)
    if exp_level == 'G':
        df.rename(columns={'gene_id': 'seq_id'}, inplace=True)
    exp_matrix = df.set_index('seq_id')
    if "all" not in geneset_id.lower():
        conn = db['sg_geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        if not geneset_records:
            bind_obj.set_error('geneset not found by query: %s', variables=(geneset_id), code="51008154")
        geneset = geneset_records['seq_list']
    else:
        geneset = list()
    if geneset:
        try:
            exp_matrix = exp_matrix.loc[geneset, :]
        except:
            bind_obj.set_error("意外错误:基因集包含基因不在该项目中")
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def get_gene_detail_whole(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip()
    annot_table = db['sg_exp']
    annot_main = annot_table.find_one({"task_id": task_id, "exp_level": "T" , 'status': 'end'})
    # if not annot_main:
    #     annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    # if not annot_main:
    #     bind_obj.set_error("Not Found in annotation_query by query %s", variables=(task_id), code="51008161")
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_exp_detail']
    query_dict = dict(exp_id=annot_main_id)
    result_dict = dict(_id=0, gene_name=1, gene_id=1, description=1, transcript_id=1)
    result = annot_detail.find(query_dict, result_dict)
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("transcript_id", inplace=True)
    #result_pd = result_pd.loc[:, ["gene_id", "gene_name", "description"]]
    #result_pd.columns = ["gene_id", "gene_name", "gene_desc"]
    result_pd = result_pd.loc[:, ["gene_id", "gene_name"]]
    result_pd.columns = ["gene_id", "gene_name"]
    result_pd["gene_name"][result_pd['gene_name'] == "-"] = None
    result_pd = result_pd.fillna(method="pad", axis=1)

    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    os.system(r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot

def export_genes_detail(data, option_name, dir_path, bind_obj=None):
    task_id, type_ = data['task_id'], data['type']
    col = db["sg_annotation_query"]
    main_id = col.find_one({"task_id": task_id, "type": type_})['main_id']

    col = db["sg_annotation_query_detail"]
    results = col.find({"query_id": ObjectId(main_id)}, {'gene_id': 1, 'gene_name': 1, 'description': 1, '_id': 0})
    # Probe Set ID    Gene Symbol     Gene Title
    out_list = []
    line_demo = '{gene_id}\t{gene_name}\t{description}\n'
    seq_id_set = set()
    out_file = os.path.join(dir_path, 'gsea.chip')
    geneset_source = bind_obj.sheet.option('geneset_source')


    with open(out_file, 'w') as out_handler:
        out_handler.write('gene_id\tgene_name\tdescription\n')
        for dic in results:
            if geneset_source == 'msigdb' and dic["gene_name"] in [None, "None", "", "-"]:
                continue
            if "description" in dic:
                if dic["description"] in [None, "None", ""]:
                    dic["description"] = "-"
            else:
                dic["description"] = "-"
            gene_id = dic['gene_id']
            if gene_id in seq_id_set:
                continue
            seq_id_set.add(gene_id)
            # bind_obj.logger.info("data is {}".format(dic))
            out_handler.write(line_demo.format(gene_id=dic["gene_id"], gene_name=dic["gene_name"].upper(), description=dic["description"]))

    return out_file

def export_genes_detail_medical(data, option_name, dir_path, bind_obj=None):
    task_id, level = data['task_id'], data['level']
    col = db["sg_exp"]
    main_id = col.find_one({"task_id": task_id, "level": level, 'is_rmbe': False})['main_id']

    col = db["sg_exp_detail"]
    if level == "G":
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
            if level == "T":
                gene_id = dic['transcript_id']
            else:
                gene_id = dic['gene_id']
            if gene_id in seq_id_set:
                continue
            seq_id_set.add(gene_id)
            # bind_obj.logger.info("data is {}".format(dic))
            if level == "G":
                try:
                    out_handler.write(line_demo.format(gene_id=dic["gene_id"], gene_name=dic["gene_name"], description=dic["description"]))
                except:
                    print dic["gene_name"]
            if level == "T":
                out_handler.write(line_demo.format(gene_id=dic["transcript_id"], gene_name=dic["gene_name"], description=dic["description"]))

    return out_file



def export_exp_matrix(data, option_name, dir_path, bind_obj=None):
    exp_id = data
    group_dict = json.loads(bind_obj.sheet.option('group'), object_pairs_hook=OrderedDict)
    geneset_id = bind_obj.sheet.option('geneset_id')
    if geneset_id.lower() in ["all", "refall"]:
        pass
    else:
        seq_list = db['sg_geneset_detail'].find_one({'geneset_id': ObjectId(geneset_id)})['seq_list']

    samples = [i for i in chain(*group_dict.values())]

    col = db["sg_exp_detail"]
    if geneset_id.lower() == "refall":
        results = col.find({"exp_id": ObjectId(exp_id), 'is_new': False}, {'_id': 0, 'exp_id': 0, 'is_new': 0})
    else:
        results = col.find({"exp_id": ObjectId(exp_id)}, {'_id': 0, 'exp_id': 0, 'is_new': 0})
    df = pd.DataFrame(list(results))
    df.set_index('seq_id', inplace=True)
    if geneset_id.lower() in ["all", "refall"]:
        df = df[samples]
    else:
        df = df[samples].loc[df.index & set(seq_list)]
    df.index.name = 'seq_id'
    out_file = os.path.join(dir_path, 'gsea.txt')
    df.to_csv(out_file, sep='\t')

    return out_file

def export_exp_matrix_medical(data, option_name, dir_path, bind_obj=None):
    # func_name = sys._getframe().f_code.co_name
    # chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    bind_obj.logger.info("step1_export_exp_matrix_medical")
    exp_id, exp_level, is_rmbe= data.split(';')
    bind_obj.logger.info("step2_export_exp_matrix_medical")
    geneset_id = bind_obj.sheet.option('geneset_id')
    if geneset_id.lower() in ["all", "refall"]:
        pass
    else:
        bind_obj.logger.info("step3_export_exp_matrix_medical")
        seq_list = db['sg_geneset_detail'].find_one({'geneset_id': ObjectId(geneset_id)})['seq_list']
        bind_obj.logger.info("step4_export_exp_matrix_medical")
    conn = db['sg_exp_detail']
    group_dict = bind_obj.sheet.option('group')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    bind_obj.logger.info("step5_export_exp_matrix_medical")
    samples = list()
    samples_raw = list()
    samples_dict = dict()
    if is_rmbe == 'false':
        for each in group_dict:
            samples += group_dict[each]
    bind_obj.logger.info("step6_export_exp_matrix_medical")
    if is_rmbe == 'true':
        for each in group_dict:
            samples_raw += group_dict[each]
            samples += ['{}_batch'.format(x) for x in group_dict[each]]
        samples_dict = dict(zip(samples, samples_raw))
    bind_obj.logger.info("step7_export_exp_matrix_medical")
    if exp_level == 'G':
        target_cols = OrderedDict(gene_id=1, _id=0)
    if exp_level == 'T':
        target_cols = OrderedDict(transcript_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    bind_obj.logger.info("step8_export_exp_matrix_medical")
    if 'kind' not in bind_obj.sheet.options() or bind_obj.sheet.options('kind') =='all':
        exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    else:
        if bind_obj.sheet.option('kind') in ["ref", 'new']:
            exp_records = conn.find({"exp_id": ObjectId(exp_id), 'kind': bind_obj.sheet.option('kind')}, target_cols)
    bind_obj.logger.info("step9_export_exp_matrix_medical")
    df = pd.DataFrame(list(exp_records))
    bind_obj.logger.info("step10_export_exp_matrix_medical")
    if exp_level == 'T':
        df.rename(columns={'transcript_id': 'seq_id'}, inplace=True)
    if exp_level == 'G':
        df.rename(columns={'gene_id': 'seq_id'}, inplace=True)
    if is_rmbe == 'true':
        df.rename(columns=samples_dict, inplace=True)
    df.set_index('seq_id', inplace=True)
    bind_obj.logger.info("step11_export_exp_matrix_medical")
    if geneset_id.lower() in ["all", "refall"]:
        df = df[samples]
    else:
        df = df[samples].loc[df.index & set(seq_list)]
    df.index.name = 'seq_id'
    out_file = os.path.join(dir_path, 'gsea.txt')
    bind_obj.logger.info("step12_export_exp_matrix_medical")
    df.to_csv(out_file, sep='\t', header=True, index=True)

    print('success to export expression matrix')
    return out_file

def export_anno_genesets(data, option_name, dir_path, bind_obj=None):
    # {'c1': data.c1, 'c2': data.c2, 'c3': data.c3, 'task_id': data.task_id,
    #                                   'type': 'origin'}
    bind_obj.logger.info("data is {}".format(data))
    task_id = bind_obj.sheet.option('task_id')
    anno_level = bind_obj.sheet.option('level')
    if (data.startswith("GO")):
        anno_type, subtype, genesets = data.split(";")
    if (data.startswith("KEGG")):
        anno_type, subtype, subtype2, genesets = data.split(";")
    if (data.startswith("DO")):
        anno_type, genesets = data.strip().split(';')
    if (data.startswith("Reactome")):
        anno_type, exp_id, genesets = data.strip().split(';')
    go = {
        'MF': 'molecular_function',
        'CC': 'cellular_component',
        'BP': 'biological_process'
    }
    kegg = {
        "Environmental Information Processing",
        "Cellular Processes",
        "Organismal Systems",
        "Human Diseases",
        "Metabolism",
        "Genetic Information Processing"
    }
    if anno_type.lower() == 'go':
        main_id = db['sg_annotation_go'].find_one({'task_id': task_id, 'type': "origin"})['main_id']
        query_dic = {'go_id': main_id, 'anno_type': anno_level}
        '''
        if subtype.upper() in go:
            query_dic['goterm'] = go[subtype.upper()]
        '''
        results = db['sg_annotation_go_list'].find(query_dic, {'gene_id': 1, 'gos_list': 1, '_id': 0})

        # if results.count() == 0:
        #    bind_obj.logger.error("go筛选结果为空, 请检查参数")
        '''
        go_list = None
        if genesets.strip():
            go_list = {i.strip() for i in genesets.strip().split(',')}
        '''
        out_file = os.path.join(dir_path, 'go.list')
        with open(out_file, 'w') as out_handler:
            out_handler.write(
                ''.join(['{}\t{}\n'.format(dic['gene_id'], dic['gos_list']) for dic in results])
            )
        return out_file
    if anno_type.lower() == 'kegg':
        main_id = db['sg_annotation_kegg'].find_one({'task_id': task_id, 'type': "origin"})['main_id']
        query_dic = {'kegg_id': main_id, 'anno_type': anno_level}
        if subtype == "defined":
            pass
        else:
            query_dic.update({'first_category': subtype})
            if subtype2 not in ['all', 'All', "","KEGG"]:
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

    if anno_type.lower() == 'reactome':
        # anno_type, genesets = data.strip().split(';')
        reactomes_list = None
        if genesets not in ["all", ""]:
            reactomes_list = {i.strip() for i in genesets.strip().split(',')}
        # main_id = db['sg_exp'].find_one({'task_id': task_id, 'level': anno_level, 'is_rmbe': False})['main_id']
        # query_dic = {'exp_id': main_id}
        query_dic = {'exp_id': ObjectId(exp_id)}
        results = db['sg_exp_detail'].find(query_dic, {'pathway_id': 1, 'pathway_description': 1, 'gene_id': 1, '_id': 0})
        print reactomes_list
        reactome_dict = dict()
        id_description = dict()
        for i in results:
            if i['pathway_id'] == "" or i['pathway_id'] == " ":
                continue
            reactome_list = list(set(i['pathway_id'].strip().split(';')))
            reactome_descriptome = list(set(i['pathway_description'].strip().split(';')))
            id_description.update(dict(zip(reactome_list, reactome_descriptome)))
            for do in reactome_list:
                if do not in reactome_dict:
                    reactome_dict[do] = i['gene_id']
                else:
                    reactome_dict[do] = reactome_dict[do] + ';' + i['gene_id']
        print reactome_dict
        out_file = os.path.join(dir_path, 'reactome.gmt')
        with open(out_file, 'w') as out_handler:
            if reactomes_list:
                out_handler.write(
                    ''.join('{}\t{}\t{}\n'.format(
                        dic, id_description[dic], reactome_dict[dic].replace(';', '\t')) for dic in reactomes_list)
                )
            else:
                out_handler.write(
                        ''.join('{}\t{}\t{}\n'.format(
                            dic, id_description[dic], reactome_dict[dic].replace(';', '\t')) for dic in reactome_dict)
                )
        return out_file

    if anno_type.lower() == 'disgenet':
        # anno_type, genesets = data.strip().split(';')
        disgenets_list = None
        if genesets not in ["all", ""]:
            disgenets_list = {i.strip() for i in genesets.strip().split(',')}
        # main_id = db['sg_exp'].find_one({'task_id': task_id, 'level': anno_level, 'is_rmbe': False})['main_id']
        # query_dic = {'exp_id': main_id}
        query_dic = {'exp_id': ObjectId(exp_id)}
        results = db['sg_exp_detail'].find(query_dic, {'pathwaw_id': 1, 'pathway_description': 1, 'gene_id': 1, '_id': 0})
        disgenet_dict = dict()
        id_description = dict()
        for i in results:
            if i['disgenet_id'] == "" or i['disgenet_id'] == " ":
                continue
            disgenet_list = list(set(i['disgenet_id'].strip().split(';')))
            disgenet_descriptome = list(set(i['disgenet_des'].strip().split(';')))
            id_description.update(dict(zip(disgenet_list, disgenet_descriptome)))
            for do in disgenet_list:
                if do not in disgenet_dict:
                    disgenet_dict[do] = i['gene_id']
                else:
                    disgenet_dict[do] = disgenet_dict[do] + ';' + i['gene_id']
        print disgenet_dict
        out_file = os.path.join(dir_path, 'disgenet.gmt')
        with open(out_file, 'w') as out_handler:
            if disgenets_list:
                out_handler.write(
                    ''.join('{}\t{}\t{}\n'.format(
                        dic, id_description[dic], disgenet_dict[dic].replace(';', '\t')) for dic in disgenets_list)
                )
            else:
                out_handler.write(
                        ''.join('{}\t{}\t{}\n'.format(
                            dic, id_description[dic], disgenet_dict[dic].replace(';', '\t')) for dic in disgenet_dict)
                )
        return out_file

    if anno_type.lower() == 'do':
        main_id = db['sg_exp'].find_one({'task_id': task_id, 'level': anno_level, "is_rmbe": False})['main_id']
        query_dic = {'exp_id': main_id}
        results = db['sg_exp_detail'].find(query_dic, {'do': 1, 'gene_id': 1, '_id': 0})
        out_file = os.path.join(dir_path, 'do.list')
        with open(out_file, 'w') as out_handler:
            out_handler.write(
                ''.join(['{}\t{}\n'.format(dic['gene_id'], ','.join([i.strip().split('(')[0] for i in dic['do'].strip().split(';')])) for dic in results])
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

    results = medical_transcriptome_db['sg_msigdb'].find(extr_dic, {'_id': 0})
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

def export_seq2name(data, option_name, dir_path, bind_obj=None):
    exp_id = data
    conn = db['sg_exp_detail']
    id2name = conn.find({"exp_id": ObjectId(exp_id), 'kind': 'ref'}, {'gene_id': 1, 'gene_name': 1, '_id': 0})
    df = pd.DataFrame(list(id2name))
    out_file = os.path.join(dir_path, 'id2name')
    df.to_csv(out_file, sep='\t', header=True, index=False)

    print('success to id2name')
    return out_file


def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} - {}'.format(k, v))
    else:
        pass
