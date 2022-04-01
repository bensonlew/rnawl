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


project_type = 'itraq_tmt'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]


def export_protein_list(data, option_name, dir_path, bind_obj=None):
    protein_list_path = os.path.join(dir_path, "%s_protein.list" % option_name)
    bind_obj.logger.debug("正在导出基因集")
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，proteinset_id:{}在sg_proteinset中未找到！".format(ObjectId(data)))
    results = collection.find_one({"proteinset_id": ObjectId(data)})
    with open(protein_list_path, "wb") as f:
        protein_list = results["seq_list"]
        for protein_id in protein_list:
            f.write(protein_id + "\n")
    return protein_list_path


def export_go_list(data, option_name, dir_path, bind_obj=None):
    go_list_path = os.path.join(dir_path, "GO.list")
    bind_obj.logger.debug("正在导出%sgo列表:%s" % (option_name, go_list_path))
    proteinset_collection = db["sg_proteinset"]
    task_id = proteinset_collection.find_one({"main_id": ObjectId(data)})["task_id"]
    my_result = db["sg_annotation_go"].find_one({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
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
            protein_id = result["gene_id"]
            go_list = result["gos_list"]
            w.write(protein_id + "\t" + go_list + "\n")
    return go_list_path

# 通过基因集id到sg_proteinset去获取task_id，type对应的信息，然后到sg_annotation_kegg去查找注释信息，然后导表；页面上的选择限制2个proteinset_id的类型
# 肯定是一样，所以任意选择一个proteinset_id过来就可以，所以接口那里随便选择了一个
# 这里的one_record是多余的，直接find的结果就可以判断，既然写了就写了吧
def export_kegg_table(data, option_name, dir_path, bind_obj=None):
    kegg_path = os.path.join(dir_path, 'protein_kegg_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg_path))
    proteinset_collection = db["sg_proteinset"]
    bind_obj.logger.debug(data)
    proteinset_result = proteinset_collection.find_one({"main_id": ObjectId(data)})
    task_id = proteinset_result["task_id"]
    bind_obj.logger.debug("ttttttt")
    bind_obj.logger.debug(task_id)
    # proteinset_type = proteinset_result["type"]
    my_result = db["sg_annotation_kegg"].find({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Protein id)\tKO_name(Protein name)\tHyperlink\tPaths\tKEGG_gene_id\n')
        for main_table in my_result:
            kegg_id = main_table["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not my_result:
                bind_obj.set_error("意外错误，annotation_kegg_id:{}在sg_annotation_kegg中未找到！".format(kegg_id))
            results = db['sg_annotation_kegg_table'].find({'kegg_id': kegg_id})
            one_record = db['sg_annotation_kegg_table'].find_one({'kegg_id': kegg_id})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:{}在sg_annotation_kegg_table中未找到！".format(ObjectId(kegg_id)))
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['ko_id'] + result['transcript_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths'], result['kegg_gene_id']))
    return kegg_path

def export_kegg_level_table(data, option_name, dir_path, bind_obj=None):
    kegg__level_path = os.path.join(dir_path, 'protein_kegg_level_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg__level_path))
    proteinset_collection = db["sg_proteinset"]
    bind_obj.logger.debug(data)
    proteinset_result = proteinset_collection.find_one({"main_id": ObjectId(data)})
    task_id = proteinset_result["task_id"]
    bind_obj.logger.debug(task_id)
    my_result = db["sg_annotation_kegg"].find({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    with open(kegg__level_path, 'wb') as w:
        w.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
        for i in my_result:
            kegg_id = i["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not kegg_id:
                bind_obj.set_error("意外错误，annotation_kegg_id:{}在sg_annotation_kegg中未找到！".format(kegg_id))
            results = db["sg_annotation_kegg_level"].find({"kegg_id": kegg_id})
            one_record = db['sg_annotation_kegg_level'].find_one({'kegg_id': kegg_id})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:{}在sg_annotation_kegg_level中未找到！".format(ObjectId(kegg_id)))
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['pathway_id'] + result['graph_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '', result['seq_num'], result['pathway_definition'],
                result['first_category'], '', result['hyperlink'], result['seq_list'], '', result['second_category']))
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
    all_list = os.path.join(dir_path, "all_protein.list")
    bind_obj.logger.debug("正在导出所有背景基因{}".format(all_list))
    collection = db['sg_express_detail']
    exp_collection = db['sg_express']
    main_collection = db['sg_proteinset']
    my_result = main_collection.find_one({"main_id": ObjectId(data)})
    task_id = my_result["task_id"]
    # proteinset_type = my_result["type"]
    bind_obj.logger.debug(task_id)
    exp_result = exp_collection.find_one({'task_id': bind_obj.sheet.option("task_id")})
    if not exp_result:
        bind_obj.set_error("意外错误，task_id:{}的背景基因在sg_proteinset中未找到！".format(data))
    exp_id = exp_result["main_id"]
    results = collection.find({"express_id": ObjectId(exp_id)})
    with open(all_list, "wb") as f:
        for result in results:
            protein_id = result['seq_id']
            f.write(protein_id + "\n")
    return all_list


def export_cog_class(data, option_name, dir_path, bind_obj=None):
    cog_path = os.path.join(dir_path, 'cog_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(cog_path))
    proteinsets, table_title, task_id, seq_type = get_proteinset_detail(data, bind_obj)
    cog_collection = db["sg_annotation_cog"]
    cog_detail_collection = db["sg_annotation_cog_detail"]
    cog_id = cog_collection.find_one({"task_id": bind_obj.sheet.option("task_id"), 'type': bind_obj.sheet.option('type')})["main_id"]
    print("cog_id:", cog_id)
    cog_results = cog_detail_collection.find({'cog_id': cog_id})
    new_table_title = []
    for tt in table_title:
        # new_tt = [tt + "_COG", tt + "_NOG", tt + "_KOG", tt + "_COG_list", tt + "_NOG_list", tt + "_KOG_list"]
        # new_tt = [tt + "_COG", tt + "_NOG", tt + "_COG_list", tt + "_NOG_list"]
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
            for gt in proteinsets:
                # kog_count = list(kog_list & proteinsets[gt][1])
                # nog_count = list(nog_list & proteinsets[gt][1])
                cog_count = list(cog_list & proteinsets[gt][1])
                # if not len(kog_count) + len(nog_count) + len(cog_count) == 0:
                #     write_line[gt] = [str(len(cog_count)), str(len(nog_count)), str(len(kog_count)), ";".join(cog_count), ";".join(nog_count), ";".join(kog_count)]
                # if not len(nog_count) + len(cog_count) == 0:
                if not len(cog_count) == 0:
                    write_line[gt] = [str(len(cog_count)), ";".join(cog_count)]
            # print('here:',  nog_list)
            print('here:',  proteinsets[gt][1])

            if len(write_line) > 0:
                w.write("{}\t{}\t".format(cr["type"], cr["function_categories"]))
                for tt in table_title:
                    w.write("\t".join(write_line[tt]) + "\t") if tt in write_line else w.write("0\tnone\t")
                w.write("\n")
    return cog_path


def get_proteinset_detail(data, bind_obj):
    proteinset_collection = db["sg_proteinset"]
    proteinsets = {}
    names = []
    task_id = ""
    proteinset_type = "P"
    for proteinset_id in data.split(","):
        proteinset_result = proteinset_collection.find_one({"main_id": ObjectId(proteinset_id)})
        if not proteinset_result:
            bind_obj.set_error("意外错误:未找到基因集_id为{}的基因集信息".format(proteinset_id))
        task_id = proteinset_result["task_id"]
        #proteinset_type = proteinset_result["type"]
        proteinset_name = proteinset_result["name"]
        names.append(proteinset_name)
        proteinsets[proteinset_name] = [proteinset_type]
        collection = db['sg_proteinset_detail']
        results = collection.find_one({"proteinset_id": ObjectId(proteinset_id)})
        proteinset_names = set(results["seq_list"])
        proteinsets[proteinset_name].append(proteinset_names)
    # print(proteinsets)
    return proteinsets, names, task_id, proteinset_type


def export_go_class(data, option_name, dir_path, bind_obj=None):
    go_path = os.path.join(dir_path, 'go_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(go_path))
    proteinsets, names, task_id, seq_type = get_proteinset_detail(data, bind_obj)
    bind_obj.logger.debug("### 正在导出{} {} {} {}".format(proteinsets, names, task_id, seq_type))
    go_collection = db["sg_annotation_go"]
    go_level_collection = db["sg_annotation_go_detail"]
    go_id = go_collection.find_one({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option('type')})["main_id"]
    one_record = go_level_collection.find_one({'go_id': go_id, "level": 2})
    if not one_record:
        bind_obj.set_error("意外错误:未找到go_id为{}的基因集信息".format(go_id))
    new_table_title = []
    for gt in proteinsets:
        new_table_title.append(gt + " num")
        new_table_title.append(gt + " percent")
        new_table_title.append(gt + " list")
    bind_obj.logger.debug(new_table_title)
    with open(go_path, "wb") as w:
        w.write("Term type\tTerm\tGO\t" + "\t".join(new_table_title) + "\n")
        term_list = ["biological_process", "cellular_component", "molecular_function"]
        for item in term_list:
            go_results = go_level_collection.find({'go_id': go_id, "level": 2})
            for gr in go_results:
                if gr["goterm"] == item:
                    seq_list = set(gr["seq_list"].split(";"))
                    write_line = {}
                    for gt in proteinsets:
                        total_protein_num = len(proteinsets[gt][1])
                        go_count = list(seq_list & proteinsets[gt][1])
                        if not len(go_count) == 0:
                            write_line[gt] = str(len(go_count)) + "\t" + str(len(go_count)/total_protein_num) + \
                                             "(" + str(len(go_count)) + "/" + str(total_protein_num) + ")" + "\t" + ";".join(go_count)
                    if len(write_line):
                        w.write("{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"]))
                        for tt in proteinsets:
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


def _get_protein_id(proteinset, proteinset_detail, _id, bind_obj=None):
    try:
        results = proteinset_detail.find_one({"proteinset_id":ObjectId(_id)})
        seq_id = results['seq_list']
    except Exception:
        bind_obj.set_error("{}在sg_proteinset_detail表中没有找到!")
    try:
        my_result = proteinset.find_one({"main_id":ObjectId(_id)})
        _name = my_result['name']
    except Exception:
        bind_obj.set_error("{}在sg_proteinset表中没有找到!")
    return seq_id, _name

# 根据"add_info":proteinset_info['task_id'] + "\t" + data.proteinset_type，就是sg_proteinset主表的task_id和页面传过来的是annotation还是annotation1
# 也就是是origin还是latest
def export_add_info(data, option_name,dir_path, bind_obj=None):
    task_id = data.split("\t")[0]
    # anno_type = data.split("\t")[1]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出add_info信息")
    col = db["sg_annotation_kegg"]
    result = col.find_one({"task_id":task_id, "type": bind_obj.sheet.option('type')})
    insert_id = result["main_id"]
    col = db["sg_annotation_kegg_level"]
    results = col.find({"kegg_id":insert_id})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info


def export_multi_protein_list(data, option_name, dir_path, bind_obj=None):
    proteinset_id = data.split(",")
    multi_proteinset_path = dir_path + "/multi_proteinset_list"
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    f = open(multi_proteinset_path, "wb")
    for n, gi in enumerate(proteinset_id):
        my_result = main_collection.find_one({'main_id': ObjectId(gi)})
        if not my_result:
            bind_obj.set_error("意外错误，proteinset_id:{}在sg_proteinset中未找到！".format(ObjectId(gi)))
        f.write(my_result["name"] + "\t")
        results = collection.find_one({"proteinset_id": ObjectId(gi)})
        f.write(",".join(results["seq_list"]) + "\n")
    return multi_proteinset_path


# ---------------------表达量相关gdq----------------------------------
def export_exp_matrix(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_express_detail']
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(accession_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    con_main = db['sg_express']
    express_id = con_main.find_one({"task_id": data, "type": "ratio"})["main_id"]

    print("999999999")
    print(target_cols)
    print(express_id) # 查找对应的id是否可以到详情表找到对应数据
    print("66666666666")
    exp_records = conn.find({"express_id": express_id}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    print("999999999")
    print(exp_matrix.head(5))
    print("66666666666")
    exp_matrix = exp_matrix.set_index('accession_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_exp_matrix2(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_express_detail']
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(accession_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    con_main = db['sg_express']
    express_id = con_main.find_one({"task_id": data, "type": "ratio"})["main_id"]

    print("999999999")
    print(target_cols)
    print(express_id) # 查找对应的id是否可以到详情表找到对应数据
    print("66666666666")
    exp_records = conn.find({"express_id": express_id}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    print("999999999")
    print(exp_matrix.head(5))
    print("66666666666")
    exp_matrix.rename(columns={'accession_id':'metab_id'}, inplace=True)
    exp_matrix = exp_matrix.set_index('metab_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_exp_matrix_scaled(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_express_detail']
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(accession_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    con_main = db['sg_express']
    express_id = con_main.find_one({"task_id": data, "type": "scaled"})["main_id"]
    exp_records = conn.find({"express_id": express_id}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('accession_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
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
        f.write('#control\tother\n')
        for each in cmp_info:
            f.write(each.replace('|', '\t') + '\n')
    return cmp_out
# ---------------------表达量相关----------------------------------


# ---------------------基因集及相关gdq----------------------------------
def export_proteinset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    exp_id, proteinset_id = data.split(";")
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
    target_cols = OrderedDict(accession_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1

    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    # get proteinset
    if not 'all' in proteinset_id:
        conn = db['sg_proteinset_detail']
        proteinset_records = conn.find_one({"proteinset_id": ObjectId(proteinset_id)})
        proteinset = proteinset_records['seq_list']
    else:
        proteinset = list()
    # get exp matrix
    conn = db['sg_express_detail']
    exp_records = conn.find({"express_id": ObjectId(exp_id)}, target_cols)

    bind_obj.logger.debug("导出表达express_id为 {}".format(exp_id))

    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('accession_id')
    if proteinset:
        filtered_proteinset = list(set(proteinset).intersection(set(exp_matrix.index)))     # added by zhangyitong on 20211108
        filtered_proteinset.sort(key=proteinset.index)
        exp_matrix = exp_matrix.loc[filtered_proteinset, :]
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_compare_exp_fc(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析accession_id 和 fc
    '''

    diff_id = bind_obj.sheet.option('diff_id')
    compare_group = bind_obj.sheet.option('compare_group')
    target_cols = OrderedDict(accession_id=1, log2fc=1, _id=0)

    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id),"compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    output = os.path.join(dir_path, option_name)
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export expression matrix')
    return output


def export_kegg_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出kegg富集表格
    '''
    kegg_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, id=1, term=1, pvalue=1, corrected_pvalue=1, seq_list=1, kegg_type=1)
    bind_obj.logger.debug("导出KEGG富集表")
    # get proteinset
    conn = db['sg_proteinset_kegg_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(kegg_enrich_id))
    kegg_enrich_records = conn.find({"kegg_enrich_id": ObjectId(kegg_enrich_id)}, target_cols)
    kegg_enrich_matrix = pd.DataFrame(list(kegg_enrich_records))
    # exp_matrix = exp_matrix.loc[proteinset, :]
    output = os.path.join(dir_path, option_name)
    kegg_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export kegg_enrich matrix')
    return output

def export_go_enrich_matrix(data, option_name, dir_path, bind_obj=None):

    go_enrich_id = bind_obj.sheet.option('enrich_id')
    go_type = bind_obj.sheet.option('go_type')

    target_cols = OrderedDict(go_id=1, go_type=1, discription=1, p_corrected=1, p_uncorrected=1, seq_list=1, depth=1, _id=0)
    bind_obj.logger.debug("导出GO富集表")
    # get proteinset
    conn = db['sg_proteinset_go_enrich_detail']
    if go_type == 'ALL' or go_type == 'All':
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id)}, target_cols)
    else:
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id), "go_type":go_type}, target_cols)
    go_enrich_matrix = pd.DataFrame(list(go_enrich_records))
    # exp_matrix = exp_matrix.loc[proteinset, :]
    output = os.path.join(dir_path, option_name)
    go_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export go_enrich matrix')
    return output

def export_protein_list_ppi(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "%s.txt" % option_name)
    bind_obj.logger.debug("正在导出蛋白集")
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，proteinset_id:{}在sg_proteinset中未找到！".format(ObjectId(data)))
    results = collection.find_one({"proteinset_id": ObjectId(data)})["seq_list"]
    with open(gene_list_path, "wb") as f:
        f.write("accession_id" + "\n")
        for result in results:
            f.write(result + "\n")
    bind_obj.logger.debug("蛋白集导出成功！")
    return gene_list_path

# for pfam stat
def export_pfam_info(data, option_name, dir_path, bind_obj=None):
    pfam_path = os.path.join(dir_path, option_name)
    need_cols = {
    "_id" : 0,
    "pfam" : 1,
    "domain" : 1,
    "description" : 1,
    "accession_id" : 1,
    "e_value" : 1,
    "length" : 1,
    "protein_start" : 1,
    "protein_end" : 1,
    "pfam_start" : 1,
    "pfam_end" : 1
}
    pfam_collection = db["sg_annotation_pfam"]
    pfam_detail_collection = db["sg_annotation_pfam_detail"]
    pfam_id = \
    pfam_collection.find_one({"task_id": bind_obj.sheet.option("task_id"), 'type': bind_obj.sheet.option('type')})[
        "main_id"]
    pfam_results = pfam_detail_collection.find({'pfam_id': pfam_id}, need_cols)
    pfam_df = pd.DataFrame(list(pfam_results))
    pfam_df.to_csv(pfam_path, sep='\t', index=False, header=True)
    return pfam_path

def export_protein_list_pfam(data, option_name, dir_path, bind_obj=None):
    protein_list_paths = list()
    proteinsets = data.split(',')
    bind_obj.logger.debug("正在导出蛋白集")
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    for proteinset in proteinsets:
        my_result = main_collection.find_one({'main_id': ObjectId(proteinset), "task_id": bind_obj.sheet.option("task_id")})
        if not my_result:
            bind_obj.set_error("意外错误，proteinset_id:{}在sg_proteinset中未找到！".format(ObjectId(proteinset)))
        name = my_result['name']
        results = collection.find_one({"proteinset_id": ObjectId(proteinset)})["seq_list"]
        protein_list_path = os.path.join(dir_path, name+'_protein.list')
        with open(protein_list_path, "wb") as f:
            for result in results:
                f.write(result + "\n")
        bind_obj.logger.debug("蛋白集%s导出成功！"%name)
        protein_list_paths.append(protein_list_path)
    return ','.join(protein_list_paths)

# for subloc stat
def export_subloc_info(data, option_name, dir_path, bind_obj=None):
    subloc_path = os.path.join(dir_path, option_name)
    need_cols = {
    "_id" : 0,
    "subloc1" : 1,
    "accession_id" : 1,
}
    subloc_collection = db["sg_annotation_subloc"]
    subloc_detail_collection = db["sg_annotation_subloc_detail"]
    subloc_id = \
    subloc_collection.find_one({"task_id": bind_obj.sheet.option("task_id"), 'type': bind_obj.sheet.option('type')})[
        "main_id"]
    subloc_results = subloc_detail_collection.find({'subloc_id': subloc_id}, need_cols)
    subloc_df = pd.DataFrame(list(subloc_results))
    subloc_df.to_csv(subloc_path, sep='\t', index=False, header=False)
    return subloc_path

# WGCNA

def download_s3_file(path, to_path):
    """
    判断文件是否在对象存储上
    """
    if os.path.exists(to_path):
        os.remove(to_path)
    elif os.path.exists(path):
        to_path = path
    elif exists(path):
        download(path, to_path)
    else:
        print('file can not find')
    return to_path

def export_wgcna_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数仅为wgcna module分析服务。
    该函数会返回wgcna_prepare得到的exp_matrix路径，
    同时还会通过查询sg_annotation_query表导出proteinid和proteinname的关系文件:seq_id2protein_name.txt
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
    query_dict = dict(query_id=annot_main_id)
    result_dict = dict(_id=0, description=1, accession_id=1)
    result = annot_detail.find(query_dict, result_dict)
    protein2name = pd.DataFrame(list(result))
    protein2name.set_index('accession_id', inplace=True)
    protein2name = protein2name.loc[list(target_seqs), :]
    protein2name.reset_index(inplace=True)
    output = os.path.join(dir_path, "seq_id2protein_name.txt")
    protein2name.to_csv(output, sep='\t', header=True, index=False)
    protein2name = pd.read_table(output, header=0)
    protein2name = protein2name.fillna(method="pad", axis=1)
    protein2name.to_csv(output, sep='\t', header=True, index=False)
    return exp_matrix+';'+output


# added  by gdq for wgcna
def export_wgcna_relate_input(data, option_name, dir_path, bind_obj=None):
    module_id = data.strip()
    # export eigenproteins
    eigenproteins = db['sg_wgcna_module_eigenproteins_detail']
    eigenproteins_found = eigenproteins.find({"module_id": ObjectId(module_id)}, {"_id": 0, "module_id":0})
    eigenproteins_pd = pd.DataFrame(list(eigenproteins_found))
    eigenproteins_pd.set_index("module", inplace=True)
    eigenproteins_path = os.path.join(dir_path, "module_eigenproteins.xls")
    eigenproteins_pd.to_csv(eigenproteins_path, sep='\t', header=True, index=True)
    # export exp matrix
    prepare_id = db['sg_wgcna_module'].find_one({"main_id": ObjectId(module_id)})['wgcna_prepare_id']
    exp_matrix = db['sg_wgcna_prepare'].find_one({"main_id": ObjectId(prepare_id)})['exp_matrix']
    # export each protein's module and protein_id/protein_name info
    membership = db['sg_wgcna_module_membership_detail']
    query_dict = dict(module_id=ObjectId(module_id))
    return_dict = dict(_id=0, seq_id=1, module=1, kme=1, block_id=1)
    result = membership.find(query_dict, return_dict)
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("seq_id", inplace=True)
    protein_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(protein_annot, sep='\t', header=True, index=True)
    return exp_matrix + ';' + eigenproteins_path + ";" + protein_annot


def export_proteinset_from_query(data, option_name, dir_path, bind_obj=None):
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
    output = os.path.join(dir_path, "all_protein.txt")
    with open(output, "w") as f:
        f.write('accession_id\n')
        for result in results:
            f.write(result['accession_id'] + '\n')
    return output

def main():
    pass

