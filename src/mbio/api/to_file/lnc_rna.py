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

project_type = 'lnc_rna'
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
        bind_obj.set_error("不能使用绝对路径或切换到其他目录!")
    if os.path.basename(to_path) == ".":
        bind_obj.set_error("目标文件不能叫\".\"!")
    target_dir = False
    if re.match(r"/$", to_path):
        target_dir = True
    work_dir = os.getcwd()
    s3transfer = S3TransferManager()
    s3transfer.base_path = work_dir
    s3transfer.overwrite = cover
    m = re.match(r"^([\w\-]+)://([\w\-]+)/(.*)$", from_file)
    if not m:
        bind_obj.set_error("下载路径%s格式不正确!" % from_file)
    else:
        region = m.group(1)
        bucket_name = m.group(2)
        key_name = m.group(3)
        if re.match(r"/$", key_name):
            if not target_dir:
                bind_obj.set_error("下载文件为文件夹时，源路径%s也必须为文件夹,以\"/\"结尾!" % to_path)
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
                        bind_obj.set_error("目标文件夹%s已经存在!" % target)
            else:
                target = os.path.join(work_dir, to_path, os.path.basename(key_name))
            s3transfer.add(from_file, target)
    s3transfer.wait()
'''

def export_gene_list(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "%s_gene.list" % option_name)
    bind_obj.logger.debug("正在导出基因集")
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(data)))
    results = collection.find_one({"geneset_id": ObjectId(data)})
    with open(gene_list_path, "wb") as f:
        gene_list = results["seq_list"]
        for gene_id in gene_list:
            f.write(gene_id + "\n")
    return gene_list_path


def export_go_list(data, option_name, dir_path, bind_obj=None):
    go_list_path = os.path.join(dir_path, "GO.list")
    bind_obj.logger.debug("正在导出%sgo列表:%s" % (option_name, go_list_path))
    geneset_collection = db["sg_geneset"]
    my_result = db["sg_annotation_go"].find_one({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    geneset_type= bind_obj.sheet.option('geneset_type')
    go_id = my_result["main_id"]
    if not my_result:
        bind_obj.set_error("意外错误，annotation_go_id:{}在sg_annotation_go中未找到！".format(go_id))
    collection = db["sg_annotation_go_list"]
    results = collection.find({"go_id": ObjectId(go_id)})
    one_record = collection.find_one({"go_id": ObjectId(go_id)})
    if not one_record:
        bind_obj.set_error("生成gos_list出错：annotation_id:{}在sg_annotation_gos_list中未找到！".format(ObjectId(go_id)))
    with open(go_list_path, "wb") as w:
        for result in results:
            gene_id = result["gene_id"]
            go_list = result["gos_list"]
            go_anno_type = result["anno_type"]
            if go_anno_type == geneset_type:
                w.write(gene_id + "\t" + go_list + "\n")
    return go_list_path

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
    my_result = db["sg_annotation_kegg"].find({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\n')
        for main_table in my_result:
            kegg_id = main_table["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not my_result:
                bind_obj.set_error("意外错误，annotation_kegg_id:{}在sg_annotation_kegg中未找到！".format(kegg_id))
            results = db['sg_annotation_kegg_table'].find({'kegg_id': kegg_id, 'anno_type': geneset_type})
            one_record = db['sg_annotation_kegg_table'].find_one({'kegg_id': kegg_id, 'anno_type': geneset_type})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:{}在sg_annotation_kegg_table中未找到！".format(ObjectId(kegg_id)))
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['ko_id'] + result['transcript_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths']))
    return kegg_path

def export_kegg_level_table(data, option_name, dir_path, bind_obj=None):
    kegg__level_path = os.path.join(dir_path, 'gene_kegg_level_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg__level_path))
    geneset_collection = db["sg_geneset"]
    bind_obj.logger.debug(data)
    geneset_result = geneset_collection.find_one({"main_id": ObjectId(data)})
    task_id = geneset_result["task_id"]
    bind_obj.logger.debug(task_id)
    my_result = db["sg_annotation_kegg"].find({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    with open(kegg__level_path, 'wb') as w:
        w.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
        for i in my_result:
            kegg_id = i["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not kegg_id:
                bind_obj.set_error("意外错误，annotation_kegg_id:{}在sg_annotation_kegg中未找到！".format(kegg_id))
            results = db["sg_annotation_kegg_level"].find({"kegg_id": kegg_id, "seq_type": "all", "anno_type": bind_obj.sheet.option("geneset_type")})
            one_record = db['sg_annotation_kegg_level'].find_one({'kegg_id': kegg_id, "seq_type": "all", 'anno_type': bind_obj.sheet.option("geneset_type")})
            if not one_record:
                results = db["sg_annotation_kegg_level"].find(
                    {"kegg_id": kegg_id, "seq_type": "ref", "anno_type": bind_obj.sheet.option("geneset_type")})
                one_record = db['sg_annotation_kegg_level'].find_one(
                    {'kegg_id': kegg_id, "seq_type": "ref", 'anno_type': bind_obj.sheet.option("geneset_type")})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:{}在sg_annotation_kegg_level中未找到！".format(ObjectId(kegg_id)))
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['pathway_id'] + result['graph_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '', result['number_of_seqs'], result['pathway_definition'],
                result['first_category'], result['anno_type'], result['hyperlink'], result['seq_list'], '', result['second_category']))
    return kegg__level_path

def export_kegg_pdf(data, option_name, dir_path, bind_obj=None):
    fs = gridfs.GridFS(db)
    annotation_collection = db["sg_annotation_kegg"]
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    main_id = annotation_collection.find_one({"task_id":task_id, "type": bind_obj.sheet.option("type")})["main_id"]
    kegg_level_collection = db["sg_annotation_kegg_level"]
    results = kegg_level_collection.find({"kegg_id":main_id, "anno_type":anno_type})
    anno_path = dir_path + "/png"
    if not os.path.exists(anno_path):
        os.mkdir(anno_path)
    for result in results:
        graph_id = result["graph_png_id"]
        pathway_id = result["pathway_id"]
        bind_obj.logger.info("正在导出{}的png文件".format(pathway_id))
        with open(anno_path + "/" + pathway_id + ".png", "w") as fw:
            content = fs.get(graph_id).read()
            fw.write(content)
    return anno_path


def export_all_list(data, option_name, dir_path, bind_obj=None):
    all_list = os.path.join(dir_path, "all_gene.list")
    bind_obj.logger.debug("正在导出所有背景基因{}".format(all_list))
    collection = db['sg_exp_detail']
    exp_collection = db['sg_exp']
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({"main_id": ObjectId(data)})
    task_id = my_result["task_id"]
    geneset_type = my_result["type"]
    bind_obj.logger.debug(task_id)
    exp_result = exp_collection.find_one({'task_id': bind_obj.sheet.option('task_id'), 'exp_level': geneset_type})
    if not exp_result:
        bind_obj.set_error("意外错误，task_id:{}的背景基因在sg_geneset中未找到！".format(data))
    exp_id = exp_result["main_id"]
    results = collection.find({"exp_id": ObjectId(exp_id)})
    with open(all_list, "wb") as f:
        for result in results:
            gene_id = result['seq_id']
            f.write(gene_id + "\n")
    return all_list


def export_cog_class(data, option_name, dir_path, bind_obj=None):
    cog_path = os.path.join(dir_path, 'cog_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(cog_path))
    genesets, table_title, task_id, geneset_type = get_geneset_detail(data, bind_obj)
    cog_collection = db["sg_annotation_cog"]
    cog_detail_collection = db["sg_annotation_cog_detail"]
    cog_id = cog_collection.find_one({"task_id": bind_obj.sheet.option('task_id'), 'type': bind_obj.sheet.option('type')})["main_id"]
    print("cog_id:", cog_id, geneset_type,  )
    if cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type}):
        cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type})
    else:
        cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "ref", 'anno_type': geneset_type})
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
                cog_count = list(cog_list & genesets[gt][1])
                # if not len(kog_count) + len(nog_count) + len(cog_count) == 0:
                #     write_line[gt] = [str(len(cog_count)), str(len(nog_count)), str(len(kog_count)), ";".join(cog_count), ";".join(nog_count), ";".join(kog_count)]
                if not len(cog_count) == 0:
                    write_line[gt] = [str(len(cog_count)), ";".join(cog_count)]
            print('here:',  genesets[gt][1])

            if len(write_line) > 0:
                w.write("{}\t{}\t".format(cr["type"], cr["function_categories"]))
                for tt in table_title:
                    w.write("\t".join(write_line[tt]) + "\t") if tt in write_line else w.write("0\tnone\t")
                w.write("\n")
    return cog_path


def get_geneset_detail(data, bind_obj):
    geneset_collection = db["sg_geneset"]
    genesets = {}
    names = []
    task_id = ""
    geneset_type = bind_obj.sheet.option('geneset_type')
    for geneset_id in data.split(","):
        geneset_result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not geneset_result:
            bind_obj.set_error("意外错误:未找到基因集_id为{}的基因集信息".format(geneset_id))
        task_id = geneset_result["task_id"]
        geneset_type = geneset_result["type"]
        geneset_name = geneset_result["name"]
        names.append(geneset_name)
        genesets[geneset_name] = [geneset_type]
        collection = db['sg_geneset_detail']
        results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
        geneset_names = set(results["seq_list"])
        genesets[geneset_name].append(geneset_names)
    # print(genesets)
    return genesets, names, task_id, geneset_type


def export_go_class(data, option_name, dir_path, bind_obj=None):
    go_path = os.path.join(dir_path, 'go_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(go_path))
    genesets, names, task_id, seq_type = get_geneset_detail(data, bind_obj)
    go_collection = db["sg_annotation_go"]
    go_level_collection = db["sg_annotation_go_detail"]
    go_id = go_collection.find_one({"task_id": bind_obj.sheet.option('task_id'), "type": bind_obj.sheet.option('type')})["main_id"]
    try:
        one_record = go_level_collection.find_one({'go_id': go_id, "seq_type": "all", "level": 2, "anno_type": seq_type})
    except:
        one_record = go_level_collection.find_one(
            {'go_id': go_id, "seq_type": "ref", "level": 2, "anno_type": seq_type})
    if not one_record:
        bind_obj.set_error("意外错误:未找到go_id为{}的基因集信息".format(go_id))
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


# ############表达量部分

def export_group_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detal这个字段
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    if data in ["all", "All", "ALL"]:
        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
        return file_path
    data = _get_objectid(data, bind_obj=bind_obj)
    group_detail = bind_obj.sheet.option('group_detail')  #另传字段
    group_table = db['sg_specimen_group']
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            bind_obj.set_error("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        bind_obj.set_error("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    group_schema = group_table.find_one({"main_id": ObjectId(data)})
    if not group_schema:
        bind_obj.set_error("无法根据传入的group_id:{}在sg_specimen_group_compare表里找到相应的记录".format(data))
    schema_name = re.sub("\s", "_", group_schema["group_name"])  # 将分组方案名的空格替换成下划线
    with open(file_path, "wb") as f:
        f.write("#sample\t" + "group" + "\n")
    sample_table_name = 'sg_specimen'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"main_id": ObjectId(sp_id)})
                if not sp:
                    bind_obj.set_error("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen_name"]
                f.write("{}\t{}\n".format(sp_name, k))
    return file_path


def _get_objectid(data, bind_obj=None):
    if not isinstance(data, ObjectId):
        if not isinstance(data, StringTypes):
            bind_obj.set_error("{}不为ObjectId类型或者其对应的字符串".format(data))
        else:
            try:
                data = ObjectId(data)
            except:
                bind_obj.set_error("{}不为ObjectId类型或者其对应的字符串".format(data))
    return data


def export_control_file(data, option_name, dir_path, bind_obj=None):  #此函数待定 不一定对
    file_path = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出计数矩阵:%s" % file_path)
    collection = db['sg_specimen_group_compare']
    result = collection.find_one({'main_id': ObjectId(data)})
    if not result:
        bind_obj.set_error("意外错误，control_id:{}在sg_specimen_group_compare中未找到！".format(ObjectId(data)))
    group_id = result['specimen_group_id']
    if group_id not in ['all', 'All', 'ALL']:
        """检查group_id的信息"""
        if isinstance(group_id, types.StringTypes):
            group_id = ObjectId(group_id)
        group_coll = db['sg_specimen_group']
        g_result = group_coll.find_one({'main_id': group_id})
        if not g_result:
            bind_obj.set_error("意外错误，control_file的group_id:{}在sg_specimen_group中未找到！".format(group_id))
    control_detail = json.loads(result['compare_names'])
    with open(file_path, 'wb') as w:
        w.write('#control\t{}\n'.format('group'))
        for i in control_detail:    #此处需要修改, 可能会有错误
            # w.write('{}\t{}\n'.format(i.keys()[0], i.values()[0]))
            control_other = i.split("|")
            w.write('{}\t{}\n'.format(control_other[0], control_other[1]))
    return file_path


def _get_gene_id(geneset, geneset_detail, _id, bind_obj=None):
    try:
        results = geneset_detail.find_one({"geneset_id":ObjectId(_id)})
        seq_id = results['seq_list']
    except Exception:
        bind_obj.set_error("{}在sg_geneset_detail表中没有找到!")
    try:
        my_result = geneset.find_one({"main_id":ObjectId(_id)})
        _name = my_result['name']
    except Exception:
        bind_obj.set_error("{}在sg_geneset表中没有找到!")
    return seq_id, _name

# 根据"add_info":geneset_info['task_id'] + "\t" + data.geneset_type，就是sg_geneset主表的task_id和页面传过来的是annotation还是annotation1
# 也就是是origin还是latest
def export_add_info(data, option_name,dir_path, bind_obj=None):
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出add_info信息")
    col = db["sg_annotation_kegg"]
    result = col.find_one({"task_id":task_id, "type": bind_obj.sheet.option('type')})
    insert_id = result["main_id"]
    col = db["sg_annotation_kegg_level"]
    results = col.find({"kegg_id":insert_id, "anno_type":anno_type})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info


def export_multi_gene_list(data, option_name, dir_path, bind_obj=None):
    geneset_id = data.split(",")
    multi_geneset_path = dir_path + "/multi_geneset_list"
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    f = open(multi_geneset_path, "wb")
    for n, gi in enumerate(geneset_id):
        my_result = main_collection.find_one({'main_id': ObjectId(gi)})
        if not my_result:
            bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(gi)))
        f.write(my_result["name"] + "\t")
        results = collection.find_one({"geneset_id": ObjectId(gi)})
        f.write(",".join(results["seq_list"]) + "\n")
    return multi_geneset_path

def export_geneset_list(data, option_name, dir_path, bind_obj=None):
    multi_geneset_path = dir_path + "/geneset_list_{}".format(data)
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']

    f = open(multi_geneset_path, "wb")
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        if data == "All":
            return "All"
        else:
            bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(gi)))
    results = collection.find_one({"geneset_id": ObjectId(data)})
    f.write("\n".join(results["seq_list"]))
    return multi_geneset_path


# ---------------------表达量相关gdq----------------------------------
def export_exp_matrix(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_exp_detail']
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    target_cols["rna_type"]=1
    target_cols["is_new"]=1
    exp_records = conn.find({"exp_id": ObjectId(data)}, target_cols)
    if 'type' in bind_obj.sheet.options():
        print("parameter 'type' found, and we will decide ref or new ")
        if 'rna_type' in bind_obj.sheet.options():
           print("parameter 'type' found, and we will decide lncRNA or mRNA or all")
           new_exp_records = list()
           if bind_obj.sheet.option('type') == "ref":
               for record in exp_records:
                   if not record['is_new']:
                       if bind_obj.sheet.option('rna_type')=="lncRNA":
                           if record['rna_type']=="lncRNA":
                             new_exp_records.append(record)
                       elif bind_obj.sheet.option('rna_type')=="mRNA":
                           if record["rna_type"]=="mRNA":
                             new_exp_records.append(record)
                       else:
                           new_exp_records.append(record)
               exp_records = new_exp_records
           elif bind_obj.sheet.option('type') == 'new':
               for record in exp_records:
                   if record['is_new']:
                       if bind_obj.sheet.option('rna_type')=="lncRNA":
                           if record['rna_type']=="lncRNA":
                             new_exp_records.append(record)
                       elif bind_obj.sheet.option('rna_type')=="mRNA":
                           if record['rna_type']=="mRNA" :
                             new_exp_records.append(record)
                       else:
                            new_exp_records.append(record)
               exp_records = new_exp_records
           else:
               for record in exp_records:
                    if bind_obj.sheet.option('rna_type') == "lncRNA":
                       if record['rna_type']=="lncRNA":
                             new_exp_records.append(record)
                    elif bind_obj.sheet.option('rna_type') == "mRNA":
                        if  record['rna_type']=="mRNA":
                          new_exp_records.append(record)
                    else:
                        new_exp_records.append(record)
               exp_records = new_exp_records
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('seq_id')
    #exp_matrix.drop(exp_matrix.columns[len(exp_matrix.columns) - 1], axis=1, inplace=True)
    exp_matrix.drop(['rna_type'], axis=1, inplace=True)
    exp_matrix.drop(['is_new'], axis=1, inplace=True)
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output


def export_count_matrix(data, option_name, dir_path, bind_obj=None):
    count_info=bind_obj.sheet.option('count_matrix')
    with open(count_info,"r") as cf:
     output = os.path.join(dir_path, option_name)
     with open(output,"w") as ff:
       if 'type' in bind_obj.sheet.options():
        print("parameter 'type' found, and we will decide ref or all ")
        if 'rna_type' in bind_obj.sheet.options():
           print("parameter 'rna_type' found, and we will decide lnc or not ")
           firstline = cf.readline()
           ff.write(firstline)
           exp_records=cf.readlines()
           if bind_obj.sheet.option('type') == "ref":
               for record in exp_records:
                   line_info = record.strip().split()
                   if not line_info[-1]=="fasle":
                       if bind_obj.sheet.option('rna_type')=="lncRNA":
                           numinfo1=len(line_info)
                           if line_info[numinfo1-2]=="lncRNA":
                               ff.write(record)
                       elif bind_obj.sheet.option('rna_type')=="mRNA":
                           numinfo1 = len(line_info)
                           if line_info[numinfo1-2] == "mRNA":
                             ff.write(record)
                       else:
                           ff.write(record)
           else:
               for record in exp_records:
                   line_info = record.strip().split()
                   if bind_obj.sheet.option('rna_type') == "lncRNA":
                       numinfo1 = len(line_info)
                       if line_info[numinfo1-2] == "lncRNA":
                           ff.write(record)
                   elif bind_obj.sheet.option('rna_type') == "mRNA":
                       numinfo1 = len(line_info)
                       if line_info[numinfo1-2] == "mRNA":
                           ff.write(record)
                   else:
                       ff.write(record)


    count_matrix = pd.read_table(output)
    count_matrix = count_matrix.set_index('seq_id')
    #exp_matrix.drop(exp_matrix.columns[len(exp_matrix.columns) - 1], axis=1, inplace=True)
    count_matrix.drop(['rna_type'], axis=1, inplace=True)
    count_matrix.drop(['is_new'], axis=1, inplace=True)
    output = os.path.join(dir_path, option_name)
    count_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export count matrix')
    return output


def export_group(data, option_name, dir_path, bind_obj=None):
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    group_out = os.path.join(dir_path, option_name)
    with open(group_out, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in group_dict:
            for each in group_dict[key]:
                f.write('{}\t{}\n'.format(each, key))
    return group_out


def export_compare(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_specimen_group_compare']
    result = conn.find_one({'main_id': ObjectId(data)})
    if not result:
        bind_obj.set_error("control_id: {} is not found in sg_specimen_group_compare".format(ObjectId(data)))
    cmp_info = json.loads(result['compare_names'])
    cmp_out = os.path.join(dir_path, option_name)
    with open(cmp_out, 'w') as f:
        f.write('#ctrl\ttest\n')
        for each in cmp_info:
            f.write(each.replace('|', '\t') + '\n')
    return cmp_out
# ---------------------表达量相关----------------------------------


# ---------------------基因集及相关gdq----------------------------------
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
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    # get geneset
    if "all" not in geneset_id.lower():
        conn = db['sg_geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        if not geneset_records:
            bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
        geneset = geneset_records['seq_list']
    else:
        geneset = list()
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
        need_indexs = exp_matrix.index & {i for i in geneset}
        exp_matrix = exp_matrix.loc[need_indexs, :]
        # exp_matrix = exp_matrix.loc[geneset, :]
    exp_matrix.dropna(how='all', inplace=True)
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_geneset_lncset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    if len(data.split(";")) == 3:
        exp_id, geneset_id, lncset_id = data.split(";")
        print "is is {} | {} | {}".format(exp_id, geneset_id, lncset_id)
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
    if geneset_id == "":
        geneset = list()
    elif "all" not in geneset_id.lower():
        conn = db['sg_geneset_detail']
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
        conn = db['sg_geneset_detail']
        lncset_records = conn.find_one({"geneset_id": ObjectId(lncset_id)})
        if not lncset_records:
            bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
        lncset = lncset_records['seq_list']
    else:
        lncset = list()
    # get all exp matrix
    conn = db['sg_exp_detail']
    if "refall" in geneset_id.lower():
        exp_records = conn.find({"exp_id": ObjectId(exp_id), 'is_new': False, 'rna_type': "mRNA"}, target_cols)
    else:
        exp_records = conn.find({"exp_id": ObjectId(exp_id), 'rna_type': "mRNA"}, target_cols)

    if "refall" in lncset_id.lower():
        lnc_exp_records = conn.find({"exp_id": ObjectId(exp_id), 'is_new': False, 'rna_type': "lncRNA"}, target_cols)
    else:
        lnc_exp_records = conn.find({"exp_id": ObjectId(exp_id), 'rna_type': "lncRNA"}, target_cols)

    exp_matrix = pd.DataFrame(list(exp_records))
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: {}'.format(exp_id))
    exp_matrix = exp_matrix.set_index('seq_id')

    lnc_exp_matrix = pd.DataFrame(list(lnc_exp_records))

    if lnc_exp_matrix.shape[0] == 0:
        # 可以没有参考lncRNA
        # bind_obj.set_error('No expression data find by query: {}'.format(exp_id))
        bind_obj.logger.info('No expression data find by query: {}'.format(exp_id))
    else:
        lnc_exp_matrix = lnc_exp_matrix.set_index('seq_id')


    print "geneset is {}".format(geneset)
    print "lncset is {}".format(lncset)

    if len(lncset) != 0 :
        lnc_exp_matrix = lnc_exp_matrix.loc[lncset, :]
    if len(geneset) != 0 :
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
    lnc_annot_table = db['sg_known_lncrna_identify']
    lncnew_annot_table = db['sg_new_lncrna_predict']
    try:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        lnc_annot_main = lnc_annot_table.find_one({"task_id": task_id, "status": "end"})
        lncnew_annot_main = lncnew_annot_table.find_one({"task_id": task_id, "status": "end"})
    except:
        bind_obj.set_error("cannot find sg_annotation_query main table")
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

    annot_detail = db['sg_annotation_query_detail']
    lnc_annot_detail = db['sg_known_lncrna_identify_detail']
    lncnew_annot_detail = db['sg_new_lncrna_predict_detail']

    exp_level = bind_obj.sheet.option('exp_level')
    if exp_level[0].upper() == 'T':
        query_dict = dict(query_id=annot_main_id,)
        result_dict = dict(_id=0, gene_name=1, gene_id=1, transcript_id=1, rna_type=1)
        result = annot_detail.find(query_dict, result_dict)
        gene2name = pd.DataFrame(list(result))
        gene2name.set_index('transcript_id', inplace=True)
        gene2name.rename(columns={"rna_type": "gene_type"}, inplace=True)
        # gene2name["type"] = "mRNA"

        '''
        lnc_query_dict = dict(query_id=annot_main_id)
        lnc_result_dict = dict(_id=0, gene_name=1, gene_id=1, lncrna_id=1)
        lnc_result = lnc_annot_detail.find(lnc_query_dict, lnc_result_dict)
        if lnc_result.count() > 0:
            lnc_gene2name = pd.DataFrame(list(lnc_result))
            lnc_gene2name.rename(columns={"lncrna_id", 'transcript_id'})
            lnc_gene2name.set_index('transcript_id', inplace=True)
            lnc_gene2name["gene_type"] = "lncRNA"
            gene2name = gene2name.append(lnc_gene2name)

        lncnew_query_dict = dict(query_id=annot_main_id)
        lncnew_result_dict = dict(id=0, gene_name=1, gene_id=1, transcript_id=1)
        lncnew_result = lncnew_annot_detail.find(lncnew_query_dict, lncnew_result_dict)
        if lncnew_result.count() > 0:
            lncnew_gene2name = pd.DataFrame(list(lncnew_result))
            lncnew_gene2name.set_index('transcript_id', inplace=True)
            lncnew_gene2name["gene_type"] = "lncRNA"
            gene2name = gene2name.append(lncnew_gene2name)
        '''

    else:
        query_dict = dict(query_id=annot_main_id, is_gene=True)
        result_dict = dict(_id=0, gene_name=1, gene_id=1, rna_type=1)
        result = annot_detail.find(query_dict, result_dict)
        gene2name = pd.DataFrame(list(result))
        gene2name.set_index('gene_id', inplace=True)
        gene2name.rename(columns={"rna_type": "gene_type"}, inplace=True)

        '''
        lnc_query_dict = dict(query_id=annot_main_id, is_gene=True)
        lnc_result_dict = dict(_id=0, gene_name=1, gene_id=1)
        lnc_result = annot_detail.find(lnc_query_dict, lnc_result_dict)
        if lnc_result.count() > 0:
            lnc_gene2name = pd.DataFrame(list(lnc_result))
            lnc_gene2name.set_index('gene_id', inplace=True)
            lnc_gene2name["gene_type"] = "lncRNA"
            gene2name = gene2name.append(lnc_gene2name)

        lncnew_query_dict = dict(query_id=annot_main_id, is_gene=True)
        lncnew_result_dict = dict(_id=0, gene_name=1, gene_id=1)
        lncnew_result = annot_detail.find(lncnew_query_dict, lncnew_result_dict)
        if lncnew_result.count() > 0:
            lncnew_gene2name = pd.DataFrame(list(lncnew_result))
            lncnew_gene2name.set_index('gene_id', inplace=True)
            lncnew_gene2name["gene_type"] = "lncRNA"
            gene2name = gene2name.append(lncnew_gene2name)
        '''

    gene2name = gene2name.loc[list(target_seqs), :]
    gene2name.reset_index(inplace=True)
    output = os.path.join(dir_path, "seq_id2gene_name.txt")
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    gene2name = pd.read_table(output, header=0)
    gene2name["gene_name"] = gene2name["gene_name"].fillna(gene2name["gene_id"])
    # gene2name.fillna(method="pad", axis=1, inplace=True)
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
            conn = db['sg_geneset_detail']
            conn_geneset = db['sg_geneset']
            geneset_records1 = conn_geneset.find_one({"main_id": ObjectId(geneset_id), "task_id": task_id})
            geneset_type = geneset_records1["type"]
            geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})

            if not geneset_records:
                bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))

            if geneset_type in geneset:
                geneset[geneset_type] = geneset[geneset_type] | set(geneset_records['seq_list'])
            else:
                geneset[geneset_type] = set(geneset_records['seq_list'])

    annot_table = db['sg_annotation_query']
    lnc_annot_table = db['sg_known_lncrna_identify']
    lncnew_annot_table = db['sg_new_lncrna_predict']

    try:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", "status": "end"})
        lnc_annot_main = lnc_annot_table.find_one({"task_id": task_id, "status": "end"})
        lncnew_annot_main = lncnew_annot_table.find_one({"task_id": task_id, "status": "end"})
    except:
        bind_obj.set_error("cannot find sg_annotation_query main table")
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

    annot_detail = db['sg_annotation_query_detail']
    lnc_annot_detail = db['sg_known_lncrna_identify_detail']
    lncnew_annot_detail = db['sg_new_lncrna_predict_detail']

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
    conn = db['sg_geneset_detail']
    conn_geneset = db['sg_geneset']
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

    genes_col = db['sg_genes']
    genes_detail_col = db['sg_genes_detail']
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
    conn = db['sg_geneset_detail']
    genes_col = db['sg_genes']
    geneset = set()
    for geneset_id in geneset_ids:
        bind_obj.logger.info("{}".format(geneset_id))
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        geneset_list = geneset_records['seq_list']
        bind_obj.logger.info("{}".format(geneset_id))
        geneset = geneset | set(geneset_list)

    # bind_obj.logger.info("gene set is {}".format(geneset))

    genes_detail_col = db['sg_genes_detail']

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
    eigengenes = db['sg_wgcna_module_eigengenes_detail']
    eigengenes_found = eigengenes.find({"module_id": ObjectId(module_id)}, {"_id": 0, "module_id":0})
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
    conn = db['sg_diff_detail']
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
    conn = db['sg_geneset_go_enrich_detail']
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
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", 'status': 'end'})
    if not annot_main:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    if not annot_main:
        task_table = db['sg_task']
        task_main = task_table.find_one({"task_id": task_id})
        if os.path.exists(task_main['annot']):
            gene_annot = os.path.join(dir_path, "seq_annot.xls")
            os.system("cut -f 1-4 " + task_main['annot'] +" |awk '{print $2\"\\t\"$1\"\\t\"$4}' > " + gene_annot)
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
    os.system(r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot


# added by gdq for tfbs predict
def export_geneid2tfid_file(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_tf_predict_detail']
    result = conn.find({"tf_predict_id": ObjectId(data.strip())}, {"_id": 0, 'gene_id': 1, "blast_hit": 1})
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("gene_id", inplace=True)
    result_file = os.path.join(dir_path, "geneid2tfid.txt")
    result_pd.to_csv(result_file, sep='\t', header=True, index=True)
    return result_file


def export_predict_result(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_tfbs_predict_detail']
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
    results = db['sg_specimen'].find({'task_id': data, 'about_qc': 'after'})
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
        rmats_info = db['sg_splicing_rmats'].find_one(
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
    geneset_detail_info = db['sg_geneset_detail'].find_one({'geneset_id': ObjectId(data)})
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
            conn = db['sg_geneset_detail']
            geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
            if not geneset_records:
                bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
            for i in geneset_records['seq_list']:
                geneset[i] = 'L' if gene_type == 'l_geneset_id' else 'G'

    output = os.path.join(dir_path, 'genes_info.json')
    with open(output, 'w') as in_handler:
        json.dump(geneset, in_handler, indent=4)

    return output

#新增批次效应的to_file函数(export_exp_batch_matrix, export_group_detail, export_exp_other_matrix)

def export_exp_batch_matrix(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    exp_id, sample_str, rna_type = data.split(';')
    sample_list = sample_str.split(',')
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in sample_list:
        target_cols[each] = 1
    collection = db['sg_exp_detail']
    print target_cols
    if rna_type.lower() == 'all':
        cursor = collection.find({'exp_id': ObjectId(exp_id)}, target_cols)
    else:
        cursor = collection.find({'exp_id': ObjectId(exp_id), 'rna_type': rna_type}, target_cols)
    output = os.path.join(dir_path, 'count.txt')
    df = pd.DataFrame(list(cursor))
    exp_matrix = df.set_index('seq_id')
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    # df.to_csv(output, sep='\t', index=False)
    return output

def export_group_detail(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    collection = db['sg_specimen_group']
    record = collection.find_one({'task_id': data})
    group = record['category_names']
    sample = record['specimen_names']
    group_dict = OrderedDict()
    for i in range(len(group)):
        group_dict[group[i]] = sample[i]
    output = os.path.join(dir_path, 'group_table.txt')
    with open(output, 'w') as f:
        f.write('#sample\tgroup\n')
        for each in group_dict:
            for key in group_dict[each]:
                f.write('{}\t{}\n'.format(key, each))
    return output

def export_exp_other_matrix(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    exp_id, exp_level = data.split(';')
    collection = db['sg_exp_detail']
    other_cols = OrderedDict(exp_id=0, _id=0)
    other = collection.find({'exp_id': ObjectId(exp_id)}, other_cols)
    output_other = os.path.join(dir_path, 'other.txt')
    df_other = pd.DataFrame(list(other))
    other_matrix = df_other.set_index('seq_id')
    other_matrix.to_csv(output_other, sep='\t', header=True, index=True)
    return output_other
