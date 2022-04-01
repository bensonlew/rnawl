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
#from biocluster.api.file.lib.s3 import S3TransferManager
#from boto.s3.bucket import Bucket
import shutil
from biocluster.api.file.lib.transfer import MultiFileTransfer
import sys

project_type = 'whole_transcriptome'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

'''
def download_from_s3(from_file, to_path="download/", cover=True):
    """
    s3 数据工作流， 来自框架workflow
    从s3对象存储下载数据到本地, 为了避免堵塞进程，此功能应该放置在流程最后执行。
    :param from_file: 需要下载的文件路径或文件路径, 必须是类似s3region://bucket/key写法。
    因为对象存储中没有文件夹的概念，需要下载文件夹必须使用"/"结尾，以明确表明下载的是文件夹
    :param to_path: 下载文件相对于当前工作目录的存放目录。
    当路径为"/"结尾时，表示下载文件存放在此文件夹下，否者为下载完整路径。
    当from_file为文件夹时，此参数也必须以"/"结尾。目录层级与下载的s3目录层级结构相同。
    默认情况下放置在当前模块工作目录的download目录下。
    :param cover: 对已存在的文件是否覆盖
    :return:
    """
    if re.match(r"^/|^\.\.", to_path):
        raise Exception("不能使用绝对路径或切换到其他目录!")
    if os.path.basename(to_path) == ".":
        raise Exception("目标文件不能叫\".\"!")
    target_dir = False
    if re.match(r"/$", to_path):
        target_dir = True
    work_dir = os.getcwd()
    s3transfer = S3TransferManager()
    s3transfer.base_path = work_dir
    s3transfer.overwrite = cover
    m = re.match(r"^([\w\-]+)://([\w\-]+)/(.*)$", from_file)
    if not m:
        raise Exception("下载路径%s格式不正确!" % from_file)
    else:
        region = m.group(1)
        bucket_name = m.group(2)
        key_name = m.group(3)
        if re.match(r"/$", key_name):
            if not target_dir:
                raise Exception("下载文件为文件夹时，源路径%s也必须为文件夹,以\"/\"结尾!" % to_path)
            conn = s3transfer.config.get_rgw_conn(region, bucket_name)
            bucket = Bucket(connection=conn, name=bucket_name)
            for key in bucket.list(prefix=key_name):
                source = os.path.join(from_file, key.name)
                target = os.path.join(target_dir, os.path.relpath(key.name, key_name))
                s3transfer.add(source, target)
        else:
            if not target_dir:  # 处理已存在文件的情况
                target = os.path.join(work_dir, to_path)
                if os.path.exists(target):
                    if cover:
                        if os.path.isdir(target):
                            shutil.rmtree(target)
                        else:
                            os.remove(target)
                    else:
                        raise Exception("目标文件夹%s已经存在!" % target)
            else:
                target = os.path.join(work_dir, to_path, os.path.basename(key_name))
            s3transfer.add(from_file, target)
    s3transfer.wait()
'''

def export_go_class(data, option_name, dir_path, bind_obj=None):
    go_path = os.path.join(dir_path, 'go_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(go_path))
    genesets, names, task_id, seq_type = get_mrna_geneset_detail(data, bind_obj)
    go_collection = db["annotation_go"]
    go_level_collection = db["annotation_go_detail"]
    go_id = go_collection.find_one({"task_id": bind_obj.sheet.option('task_id')})["main_id"]
    one_record = go_level_collection.find_one({'go_id': go_id, "level": 2, "anno_type": seq_type})
    if not one_record:
        bind_obj.set_error("意外错误:未找到go_id为%s的基因集信息", variables=(go_id), code="51008141")
    new_table_title = []
    for gt in genesets:
        new_table_title.append(gt + " num")
        new_table_title.append(gt + " percent")
        new_table_title.append(gt + " list")
    bind_obj.logger.debug(new_table_title)
    with open(go_path, "wb") as w:
        w.write("Term type\tTerm\tGO\t" + "\t".join(new_table_title) + "\n")
        term_list = ["biological_process", "cellular_component", "molecular_function"]
        for item in term_list:
            if go_level_collection.find_one({'go_id': go_id, "seq_type": "all", "level": 2, "anno_type": seq_type}) :
                go_results = go_level_collection.find({'go_id': go_id, "seq_type": "all", "level": 2, "anno_type": seq_type})
            else:
                go_results = go_level_collection.find(
                    {'go_id': go_id, "seq_type": "ref", "level": 2, "anno_type": seq_type})
            for gr in go_results:
                if gr["goterm"] == item:
                    seq_list = set(gr["seq_list"].split(";"))
                    write_line = {}
                    for gt in genesets:
                        total_gene_num = len(genesets[gt][1])
                        go_count = list(seq_list & genesets[gt][1])
                        if not len(go_count) == 0:
                            write_line[gt] = str(len(go_count)) + "\t" + str(len(go_count)/total_gene_num) + \
                                             "(" + str(len(go_count)) + "/" + str(total_gene_num) + ")" + "\t" + ";".join(go_count)
                    if len(write_line):
                        w.write("{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"]))
                        for tt in genesets:
                            w.write(write_line[tt] + "\t") if tt in write_line else w.write("0\t0\tnone\t")
                        w.write("\n")
    return go_path


def export_cog_class(data, option_name, dir_path, bind_obj=None):
    cog_path = os.path.join(dir_path, 'cog_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(cog_path))
    genesets, table_title, task_id, geneset_type = get_mrna_geneset_detail(data, bind_obj)
    bind_obj.logger.debug("正在导出{}".format("nonono"))
    cog_collection = db["annotation_cog"]
    cog_detail_collection = db["annotation_cog_detail"]
    cog_id = cog_collection.find_one({"task_id": bind_obj.sheet.option('task_id')})["main_id"]
    print("cog_id:", cog_id, geneset_type,  )
    cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type})
    # 更改mongo查询是否为空的判断
    num = 0
    for i in cog_results:
        num += 1
    if num == 0:
        cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "ref", 'anno_type': geneset_type})
    else:
        # cog_results循环完了需要重新赋值，否则为空
        cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type})
    # if cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type}):
    #     cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type})
    # else:
    #     cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "ref", 'anno_type': geneset_type})
    new_table_title = []
    for tt in table_title:
        # new_tt = [tt + "_COG", tt + "_NOG", tt + "_KOG", tt + "_COG_list", tt + "_NOG_list", tt + "_KOG_list"]
        new_tt = [tt + "_COG", tt + "_COG_list"]
        new_table_title = new_table_title + new_tt
    bind_obj.logger.debug(table_title)
    with open(cog_path, "wb") as w:
        w.write("Type\tFunctional Categoris\t" + "\t".join(new_table_title) + "\n")
        for cr in cog_results:
            # kog_list = set([])
            # nog_list = set(cr["nog_list"].split(";") if cr["nog_list"] else [])
            cog_list = set(cr["cog_list"].split(";") if cr["cog_list"] else [])
            write_line = {}
            for gt in genesets:
                # kog_count = list(kog_list & genesets[gt][1])
                # nog_count = list(nog_list & genesets[gt][1])
                print gt
                print genesets[gt][1]
                print cog_list
                cog_count = list(cog_list & genesets[gt][1])
                # if not len(kog_count) + len(nog_count) + len(cog_count) == 0:
                #     write_line[gt] = [str(len(cog_count)), str(len(nog_count)), str(len(kog_count)), ";".join(cog_count), ";".join(nog_count), ";".join(kog_count)]
                if not len(cog_count) == 0:
                    write_line[gt] = [str(len(cog_count)), ";".join(cog_count)]
            print('here:',  genesets[gt][1])

            if len(write_line) > 0:
                w.write("{}\t{}\t".format(cr["type"], cr["function_categories"]))
                # all_line = ""
                for tt in table_title:
                    # line = "\t".join(write_line[tt]) + "\t" if  tt in write_line else "0\tnone\t"
                    # all_line += line
                    w.write("\t".join(write_line[tt]) + "\t") if tt in write_line else w.write("0\tnone\t")
                w.write("\n")
                # all_line = all_line.strip()
                # w.write(all_line + "\n")
    return cog_path

def get_geneset_detail(data, bind_obj):
    geneset_collection = db["geneset"]
    genesets = {}
    names = []
    task_id = ""
    geneset_type = bind_obj.sheet.option('geneset_type')
    for geneset_id in data.split(","):
        geneset_result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not geneset_result:
            bind_obj.set_error("意外错误:未找到基因集_id为%s的基因集信息", variables=(geneset_id), code="51008140")
        task_id = geneset_result["task_id"]
        geneset_type = geneset_result["type"]
        geneset_name = geneset_result["name"]
        names.append(geneset_name)
        genesets[geneset_name] = [geneset_type]
        collection = db['geneset_detail']
        results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
        geneset_names = set(results["seq_list"])
        genesets[geneset_name].append(geneset_names)
    #print genesets
    return genesets, names, task_id, geneset_type

def get_mrna_geneset_detail(data, bind_obj):
    geneset_collection = db["geneset"]
    genesets = {}
    names = []
    task_id = ""
    geneset_type = bind_obj.sheet.option('geneset_type')
    for geneset_id in data.split(","):
        geneset_result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not geneset_result:
            bind_obj.set_error("意外错误:未找到基因集_id为%s的基因集信息", variables=(geneset_id), code="51008140")
        task_id = geneset_result["task_id"]
        geneset_type = geneset_result["level"]
        geneset_name = geneset_result["name"]
        names.append(geneset_name)
        genesets[geneset_name] = [geneset_type]
        collection = db['geneset_detail']
        results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
        if 'category_list' not in results:
            bind_obj.set_error(
                'geneset_detail table must contain "category_list" field when provide one geneset id')
        seq_list = results["seq_list"]
        category_list = results["category_list"]
        mrna_list = []
        for seq_id, category in zip(seq_list, category_list):
            if category == "mRNA":
                mrna_list.append(seq_id)
        geneset_names = set(mrna_list)
        genesets[geneset_name].append(geneset_names)
    # print genesets
    return genesets, names, task_id, geneset_type

def export_multi_gene_list(data, option_name, dir_path, bind_obj=None):
    data = json.loads(data)
    geneset_id_list = data['geneset_id'].split(",")
    source = data['source']
    multi_geneset_path = dir_path + "/multi_geneset_list"
    main_collection = db['geneset']
    collection = db['geneset_detail']
    with open(multi_geneset_path, "wb") as out_handler:
        if len(geneset_id_list) == 1 and source == 'diff':
            geneset_id = geneset_id_list[0]
            my_result = main_collection.find_one({'main_id': ObjectId(geneset_id)})
            if not my_result:
                bind_obj.set_error("意外错误，geneset_id:%s在geneset中未找到！", variables=(ObjectId(geneset_id)), code="55700108")
            geneset_name = my_result["name"]
            results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
            seq_list = results["seq_list"]
            if 'regulate_list' not in results:
                bind_obj.set_error(
                    'geneset_detail table must contain "regulate_list" field when provide one geneset id')
            if 'category_list' not in results:
                bind_obj.set_error(
                    'geneset_detail table must contain "category_list" field when provide one geneset id')
            regulate_list = results["regulate_list"]
            category_list= results["category_list"]
            up_mrna_list, down_mrna_list = [], []
            for seq_id, regulate,category in zip(seq_list, regulate_list,category_list):
                if regulate == 'up' and category == "mRNA":
                    up_mrna_list.append(seq_id)
                elif category == "mRNA" :
                    down_mrna_list.append(seq_id)
            out_handler.write(geneset_name + '_up\t' + ','.join(up_mrna_list) + '\n')
            out_handler.write(geneset_name + '_down\t' + ','.join(down_mrna_list) + '\n')
        else:
            for n, gi in enumerate(geneset_id_list):
                my_result = main_collection.find_one({'main_id': ObjectId(gi)})
                if not my_result:
                    bind_obj.set_error("意外错误，geneset_id:%s在geneset中未找到！", variables=(ObjectId(gi)), code="55700110")
                out_handler.write(my_result["name"] + "\t")
                results = collection.find_one({"geneset_id": ObjectId(gi)})
                out_handler.write(",".join(results["seq_list"]) + "\n")

    return multi_geneset_path

# 通过基因集id到geneset去获取task_id，type对应的信息，然后到annotation_kegg去查找注释信息，然后导表；页面上的选择限制2个geneset_id的类型
# 肯定是一样，所以任意选择一个geneset_id过来就可以，所以接口那里随便选择了一个
# 这里的one_record是多余的，直接find的结果就可以判断，既然写了就写了吧
def export_kegg_table(data, option_name, dir_path, bind_obj=None):
    kegg_path = os.path.join(dir_path, 'gene_kegg_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg_path))
    geneset_collection = db["geneset"]
    bind_obj.logger.debug(data)
    geneset_result = geneset_collection.find_one({"main_id": ObjectId(data)})
    task_id = geneset_result["task_id"]
    bind_obj.logger.debug("ttttttt")
    bind_obj.logger.debug(task_id)
    geneset_type = geneset_result["level"]
    my_result = db["annotation_kegg"].find({"task_id": bind_obj.sheet.option("task_id")})
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\n')
        for main_table in my_result:
            kegg_id = main_table["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not my_result:
                bind_obj.set_error("意外错误，annotation_kegg_id:%s在annotation_kegg中未找到！", variables=(kegg_id), code="51008135")
            results = db['annotation_kegg_table'].find({'kegg_id': kegg_id, 'anno_type': geneset_type})
            one_record = db['annotation_kegg_table'].find_one({'kegg_id': kegg_id, 'anno_type': geneset_type})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:%s在annotation_kegg_table中未找到！", variables=(ObjectId(kegg_id)), code="51008136")
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['ko_id'] + result['transcript_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths']))
    return kegg_path

def export_kegg_level_table(data, option_name, dir_path, bind_obj=None):
    kegg__level_path = os.path.join(dir_path, 'gene_kegg_level_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg__level_path))
    geneset_collection = db["geneset"]
    bind_obj.logger.debug(data)
    geneset_result = geneset_collection.find_one({"main_id": ObjectId(data)})
    task_id = geneset_result["task_id"]
    bind_obj.logger.debug(task_id)
    my_result = db["annotation_kegg"].find({"task_id": bind_obj.sheet.option("task_id")})
    with open(kegg__level_path, 'wb') as w:
        w.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
        for i in my_result:
            kegg_id = i["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not kegg_id:
                bind_obj.set_error("意外错误，annotation_kegg_id:%s在annotation_kegg中未找到！", variables=(kegg_id), code="51008137")
            results = db["annotation_kegg_level"].find({"kegg_id": kegg_id, "seq_type": "all", "anno_type": bind_obj.sheet.option("geneset_type")})
            one_record = db['annotation_kegg_level'].find_one({'kegg_id': kegg_id, "seq_type": "all", 'anno_type': bind_obj.sheet.option("geneset_type")})
            if not one_record:
                results = db["annotation_kegg_level"].find(
                    {"kegg_id": kegg_id, "seq_type": "ref", "anno_type": bind_obj.sheet.option("geneset_type")})
                one_record = db['annotation_kegg_level'].find_one(
                    {'kegg_id': kegg_id, "seq_type": "ref", 'anno_type': bind_obj.sheet.option("geneset_type")})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:%s在annotation_kegg_level中未找到！", variables=(ObjectId(kegg_id)), code="51008138")
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['pathway_id'] + result['graph_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '', result['number_of_seqs'], result['pathway_definition'],
                result['first_category'], result['anno_type'], result['hyperlink'], result['seq_list'], '', result['second_category']))
    return kegg__level_path

# 根据"add_info":geneset_info['task_id'] + "\t" + data.geneset_type，就是geneset主表的task_id和页面传过来的是annotation还是annotation1
# 也就是是origin还是latest
def export_add_info(data, option_name,dir_path, bind_obj=None):
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出add_info信息")
    col = db["annotation_kegg"]
    result = col.find_one({"task_id":task_id})
    insert_id = result["main_id"]
    col = db["annotation_kegg_level"]
    results = col.find({"kegg_id":insert_id, "anno_type":anno_type})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info


def export_gene_list(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "%s_gene.list" % option_name)
    bind_obj.logger.debug("正在导出基因集")
    collection = db['geneset_detail']
    main_collection = db['geneset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:%s在geneset中未找到！", variables=(ObjectId(data)), code="51008132")
    results = collection.find_one({"geneset_id": ObjectId(data)})
    with open(gene_list_path, "wb") as f:
        gene_list = results["seq_list"]
        for gene_id in gene_list:
            f.write(gene_id + "\n")
    return gene_list_path


def export_mrna_gene_list(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "%s_gene.list" % option_name)
    bind_obj.logger.debug("正在导出基因集")
    collection = db['geneset_detail']
    main_collection = db['geneset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:%s在geneset中未找到！", variables=(ObjectId(data)), code="51008132")
    results = collection.find_one({"geneset_id": ObjectId(data)})
    with open(gene_list_path, "wb") as f:
        gene_list = results["seq_list"]

        for gene_id in gene_list:
            f.write(gene_id + "\n")
    return gene_list_path


def export_all_list(data, option_name, dir_path, bind_obj=None):
    all_list = os.path.join(dir_path, "all_gene.list")
    bind_obj.logger.debug("正在导出所有背景基因{}".format(all_list))
    collection = db['exp_detail']
    exp_collection = db['exp']
    main_collection = db['geneset']
    my_result = main_collection.find_one({"main_id": ObjectId(data)})
    task_id = my_result["task_id"]
    geneset_type = my_result["level"]
    bind_obj.logger.debug(task_id)
    exp_result = exp_collection.find_one({'task_id': bind_obj.sheet.option('task_id'), 'level': geneset_type})
    if not exp_result:
        bind_obj.set_error("意外错误，task_id:%s的背景基因在geneset中未找到！", variables=(data), code="51008139")
    exp_id = exp_result["main_id"]
    # results = collection.find({"exp_id": ObjectId(exp_id)})
    #测试数据有部分数据不对应，缺少字段，因此临时用该代码，上线后需取消transcript_id:{$exists:false}
    results = collection.find({"exp_id": ObjectId(exp_id),"transcript_id":{"$exists":True}})
    if geneset_type == "G":
        seq_type = "gene_id"
    else:
        seq_type = "transcript_id"
    with open(all_list, "wb") as f:
        for result in results:
            gene_id = result[seq_type]
            f.write(gene_id + "\n")
    return all_list

def export_go_list(data, option_name, dir_path, bind_obj=None):
    go_list_path = os.path.join(dir_path, "GO.list")
    bind_obj.logger.debug("正在导出%sgo列表:%s" % (option_name, go_list_path))
    geneset_collection = db["geneset"]
    # task_id = geneset_collection.find_one({"main_id": ObjectId(data)})["task_id"]
    my_result = db["annotation_go"].find_one({"task_id": bind_obj.sheet.option("task_id")})
    geneset_type= bind_obj.sheet.option('geneset_type')
    go_id = my_result["main_id"]
    if not my_result:
        bind_obj.set_error("意外错误，annotation_go_id:%s在annotation_go中未找到！", variables=(go_id), code="51008133")
    collection = db["annotation_go_list"]
    results = collection.find({"go_id": ObjectId(go_id)})
    one_record = collection.find_one({"go_id": ObjectId(go_id)})
    if not one_record:
        bind_obj.set_error("生成gos_list出错：annotation_id:%s在annotation_gos_list中未找到！", variables=(ObjectId(go_id)), code="51008134")
    with open(go_list_path, "wb") as w:
        for result in results:
            gene_id = result["gene_id"]
            go_list = result["gos_list"]
            go_anno_type = result["anno_type"]
            if go_anno_type == geneset_type:
                w.write(gene_id + "\t" + go_list + "\n")
    return go_list_path


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

def get_gene_detail(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip()
    annot_table = db['exp']
    annot_main = annot_table.find_one({"task_id": task_id, "level": "T" , 'status': 'end'})
    # if not annot_main:
    #     annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    # if not annot_main:
    #     bind_obj.set_error("Not Found in annotation_query by query %s", variables=(task_id), code="51008161")
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['exp_detail']
    query_dict = dict(exp_id=annot_main_id,category="mRNA")
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

def get_gene_detail_whole(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip()
    annot_table = db['exp']
    annot_main = annot_table.find_one({"task_id": task_id, "level": "T" , 'status': 'end'})
    # if not annot_main:
    #     annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    # if not annot_main:
    #     bind_obj.set_error("Not Found in annotation_query by query %s", variables=(task_id), code="51008161")
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['exp_detail']
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

def get_gene_type(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip().split(",")[0]
    level = data.strip().split(",")[1]
    category_table = db['exp']
    category_main = category_table.find_one({"task_id": task_id, "level": level , 'status': 'end'})
    # if not annot_main:
    #     annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    # if not annot_main:
    #     bind_obj.set_error("Not Found in annotation_query by query %s", variables=(task_id), code="51008161")
    if "main_id" not in category_main:
        category_main_id = category_main['_id']
    else:
        category_main_id = category_main['main_id']
    category_detail = db['exp_detail']
    query_dict = dict(exp_id=category_main_id, )
    if level == "G":
        result_dict = dict(_id=0, gene_id=1, category=1)
    else:
        result_dict = dict(_id=0, transcript_id=1, category=1)
    result = category_detail.find(query_dict, result_dict)
    result_pd = pd.DataFrame(list(result))
    if level == "G":
        result_pd.set_index("gene_id", inplace=True)
    else:
        result_pd.set_index("transcript_id", inplace=True)
    gene_category = os.path.join(dir_path, "seq_category.xls")
    result_pd.to_csv(gene_category, sep='\t', header=True, index=True)
    # os.system(r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_category


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

def export_geneset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    exp_id, geneset_id = data.split(";")
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
    geneset_type=bind_obj.sheet.option('gt')
    if geneset_type == "G":
        seq_type = "gene_id"
        target_cols = OrderedDict(gene_id=1, _id=0)
    else:
        seq_type = "transcript_id"
        target_cols = OrderedDict(transcript_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    # get geneset
    geneset = list()
    if "all" not in geneset_id.lower():
        conn = db['geneset_detail']
        for i in geneset_id.split(","):
            geneset_records = conn.find_one({"geneset_id": ObjectId(i)})
            if not geneset_records:
                bind_obj.set_error('geneset not found by query: %s', variables=(geneset_id), code="51008154")
            genesetp = geneset_records['seq_list']
            geneset.extend(genesetp)
    else:
        geneset = list()
    geneset = list(set(geneset))
    # get all exp matrix
    conn = db['exp_detail']
    if "refall" in geneset_id.lower():
        exp_records = conn.find({"exp_id": ObjectId(exp_id), 'is_new': False}, target_cols)
    else:
        exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: %s', variables=(exp_id), code="51008155")
    exp_matrix=exp_matrix.rename(columns={seq_type:"seq_id"})
    exp_matrix = exp_matrix.set_index('seq_id')
    if geneset:
        exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output


def export_kegg_enrich_info_filter(data,option_name,dir_path,bind_obj=None):
    collection = db["geneset_kegg_enrich_detail"]
    target_cols = OrderedDict(id=1, _id=0, first_category=1, second_category=1, seq_list=1)
    new_enrich_records = list()
    if bind_obj.sheet.option("stat_level"):
        print("parameter 'stat_level' found, and we will filter data by {}").format(bind_obj.sheet.option("stat_level"))
    if bind_obj.sheet.option("stat_threshold_value") and bind_obj.sheet.option("stat_numbers_value"):
        print("stat_threshold_value is {},stat_numbers_value is {}").format(bind_obj.sheet.option("stat_threshold_value"),bind_obj.sheet.option("stat_numbers_value"))
    v1=float(bind_obj.sheet.option("stat_threshold_value"))
    v2=int(bind_obj.sheet.option("stat_numbers_value"))
    if bind_obj.sheet.option("stat_level")=="pvalue":
        target_cols["pvalue"] = 1
        enrich_records = collection.find({"kegg_enrich_id": ObjectId(data)}, target_cols)
        for record in enrich_records:
           if  float(record["pvalue"]) <= v1:
               new_enrich_records.append(record)
        enrich_records=new_enrich_records
    else:
        target_cols["corrected_pvalue"] = 1
        enrich_records = collection.find({"kegg_enrich_id": ObjectId(data)}, target_cols)
        for record in enrich_records:
            if record["corrected_pvalue"] <= v1:
                new_enrich_records.append(record)
        enrich_records = new_enrich_records
    if len(enrich_records) < v2 :
        bind_obj.set_error("filted_records less than %s ", variables=(bind_obj.sheet.option("stat_numbers_value")), code="51008157")
    enrich_matrix = pd.DataFrame(list(enrich_records))
    enrich_matrix = enrich_matrix.set_index('id')
    enrich_matrix=enrich_matrix.sort_values(bind_obj.sheet.option("stat_level"),ascending=True)[:v2]
    output = os.path.join(dir_path, option_name)
    enrich_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export geneset enrich filterd matrix')
    return output


def export_kegg_enrich_info_ids(data,option_name,dir_path,bind_obj=None):
    collection=db["sg_geneset_kegg_enrich_detail"]
    if 'pathway_ids' in bind_obj.sheet.options():
        print("parameter 'pathway_ids' found, and we will choose ids for id_list ")
    ids=bind_obj.sheet.option('pathway_ids').split(",")
    target_cols = OrderedDict(id=1, _id=0,first_category=1,second_category=1,seq_list=1)
    enrich_records = collection.find({"kegg_enrich_id": ObjectId(data)}, target_cols)
    new_enrich_records = list()
    for record in enrich_records:
        if record['id'] in ids :
            new_enrich_records.append(record)
    enrich_records = new_enrich_records
    if len(enrich_records) == 0:
        bind_obj.set_error("filted_records is empty", code="51008156")
    enrich_matrix = pd.DataFrame(list(enrich_records))
    enrich_matrix = enrich_matrix.set_index('id')
    output = os.path.join(dir_path, option_name)
    enrich_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export geneset kegg enrich detail matrix')
    return output

