# -*- coding: utf-8 -*-

import os
import re
import json
from types import StringTypes
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import defaultdict
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.controllers.project.metabolome_controller import MetabolomeController


client = Config().get_mongo_client(mtype="metabolome")
db = client[Config().get_mongo_dbname("metabolome")]

def _get_objectid(data):
    if not isinstance(data, ObjectId):
        if not isinstance(data, StringTypes):
            raise Exception("{}不为ObjectId类型或者其对应的字符串".format(data))
        else:
            try:
                data = ObjectId(data)
            except:
                raise Exception("{}不为ObjectId类型或者其对应的字符串".format(data))
    return data

def export_metab_set(data, option_name, dir_path, bind_obj=None):   # add by haidong.gu
    """
    根据代谢集主表_id获取代谢集列表
    """
    file_path = os.path.join(dir_path, "%s_input.set.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的代谢集为文件，路径：%s" % (option_name, file_path))
    data = _get_objectid(data)
    main_id = db['metab_set'].find_one({"_id": data})["main_id"]
    main_id = _get_objectid(main_id)
    metab_set_list = db['metab_set_detail'].find_one({"set_id": main_id})['set_list']
    #metab_set_list = db['metab_set_detail'].find_one({"set_id": data})['set_list']
    out_file = open(file_path,"w")
    for i in metab_set_list:
        out_file.write(i + "\n")
    out_file.close()
    return file_path

def export_mul_metab_set(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s_input.set.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的代谢集为文件，路径：%s" % (option_name, file_path))
    f = open(file_path, "w")
    for i in data.split(","):
        set_id = _get_objectid(i)
        main_id = db['metab_set'].find_one({"_id": set_id})["main_id"]
        main_id = _get_objectid(main_id)
        metab_set_list = db['metab_set_detail'].find_one({"set_id": main_id})['set_list']
        #metab_set_list = db["metab_set_detail"].find_one({"set_id": set_id})["set_list"]
        name = db["metab_set"].find_one({"_id": set_id})["name"]
        f.write(name + "\t" + ",".join(metab_set_list) + "\n")
    f.close()
    return  file_path

def export_metab_detail(data, option_name, dir_path, bind_obj=None):
    # 直接用本地文件，不导表
    pass

def export_overview(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s_input.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s为文件，路径：%s" % (option_name, file_path))
    f = open(file_path, 'w')
    main_id = db['anno_overview'].find_one({"task_id": data})["main_id"]
    # data = _get_objectid(data)
    result = db['anno_overview_detail'].find({"overview_id": main_id})
    head_list = ["metab_id", "metab", "mode", "hmdb_id",  "cas_id", "formula",
                      "compound_id", "compound_name", "compound_first_category", "compound_second_category",
                      "pathway_id", "description", "kegg_first_category", "kegg_second_category"]  #rm element hmdb_name "mass",
    f.write("\t".join(head_list) + "\n")
    for one in result:
        mongo_name_list = ['metab_id', 'metab', 'mode', 'hmdb_id',  'cas_id', 'formula',
                           'compound_id', 'compound_name', 'c1_category', 'c2_category', 'pathway_id', 'description',
                           'p1_category', 'p2_category']  #rm element  hmdb_name 'mass',
        # c_mongo_name_list = []
        # for c in mongo_name_list:
        #     if c in one.keys():
        #         c_mongo_name_list.append(c)
        write_list = map(lambda x: one[x], mongo_name_list)
        f.write("\t".join(write_list) + "\n")
    f.close()
    return file_path

def export_overview_ko(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s_input.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s为文件， 路径:%s" % (option_name, file_path))
    f = open(file_path, 'w')
    main_id = db['anno_overview'].find_one({"task_id": data})["main_id"]
    # data = _get_objectid(data)
    result = db['anno_overview_ko'].find({"overview_id": main_id})
    head_list = ["pathway_id", "description", "first_category", "second_category", "compound_id", "metab_id", "hyperlink"]
    f.write("\t".join(head_list) + "\n")
    for one in result:
        write_list = map(lambda x: one[x], head_list)
        f.write("\t".join(write_list) + "\n")
    f.close()
    return file_path

def export_group_by_detail(data, option_name, dir_path, bind_obj=None):  # add by shaohua.yuan
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detail这个字段
    data 为group_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    #if data not in ["all","All","ALL"]:
    #data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    if not isinstance(group_detail, dict):
        try:
            bind_obj.logger.debug("正确")
            table_dict = json.loads(group_detail)
            bind_obj.logger.debug(table_dict)
        except Exception:
            bind_obj.logger.debug("错误")
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    else:
        table_dict = group_detail
    if not isinstance(table_dict, dict):
        bind_obj.logger.debug("错误")
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    schema_name = "group_name"
    with open(file_path, "wb") as f:
        f.write("#sample\t" + schema_name + "\n")
        for k in table_dict:
            sample_names = table_dict[k]
            for each in sample_names:
                if len(table_dict) == 1 and k in ['all', 'All', 'ALL']:
                    f.write("{}\t{}\n".format(each, each))
                else:
                    f.write("{}\t{}\n".format(each, k))
    return file_path

def export_group_by_detail2(data, option_name, dir_path, bind_obj=None):  # add by shaohua.yuan
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detail这个字段
    data 为group_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    #if data not in ["all","All","ALL"]:
    #data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    if not isinstance(group_detail, dict):
        try:
            bind_obj.logger.debug("正确")
            table_dict = json.loads(group_detail)
            bind_obj.logger.debug(table_dict)
        except Exception:
            bind_obj.logger.debug("错误")
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    else:
        table_dict = group_detail
    if not isinstance(table_dict, dict):
        bind_obj.logger.debug("错误")
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    schema_name = "group_name"
    with open(file_path, "wb") as f:
        f.write("#sample\t" + schema_name + "\n")
        for k in table_dict:
            sample_names = table_dict[k]
            for each in sample_names:
                if len(table_dict) == 1 and k in ['all', 'All', 'ALL']:
                    f.write("{}\t{}\n".format(each, k))  # export_group_by_detail2 和export_group_by_detail的不同
                else:
                    f.write("{}\t{}\n".format(each, k))
    return file_path


def export_group_table(data, option_name, dir_path, bind_obj=None):  # add by shaohua.yuan
    """
    按group_id获取带QC分组的group文件
    使用时确保你的workflow的option里group_detail这个字段
    data 为task_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件
    """
    file_path = os.path.join(dir_path, "%s_input.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data = _get_objectid(data)
    # group_detail = bind_obj.sheet.option('group_detail')
    group_info = db["specimen_group"].find_one({"_id":ObjectId(data)})
    category_names = group_info["category_names"]
    specimen_names = group_info["specimen_names"]
    metabolome = Metabolome()
    #id_to_sample = id2name(task_id, type="task")
    with open(file_path, "wb") as f:
        f.write("#sample\tgroup_name\n")
        for i in range(0,len(category_names)):
            group_specimen = specimen_names[i]
            for each in group_specimen:
                f.write("{}\t{}\n".format(each, category_names[i]))
    return file_path

def export_group_table_noQC(data, option_name, dir_path, bind_obj=None):  # add by shaohua.yuan
    """
    按group_id获取不带QC分组的group文件
    使用时确保你的workflow的option里group_detail这个字段
    data 为task_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件
    """
    file_path = os.path.join(dir_path, "%s_input.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data = _get_objectid(data)
    #group_detail = bind_obj.sheet.option('group_detail')
    group_info = db["specimen_group"].find_one({"_id":ObjectId(data)})
    category_names = group_info["category_names"]
    specimen_names = group_info["specimen_names"]
    metabolome = Metabolome()
    #id_to_sample = id2name(task_id, type="task")
    with open(file_path, "wb") as f:
        f.write("#sample\tgroup_name\n")
    with open(file_path, "ab") as f:
        for i in range(0,len(category_names)):
            if category_names[1] != "QC":
                group_specimen = specimen_names[i]
                for each in group_specimen:
                    f.write("{}\t{}\n".format(each, category_names[i]))
    return file_path

def export_metab_set1(data, option_name, dir_path, bind_obj=None):
    """
    根据代谢集主表_id获取代谢集列表
    """
    file_path = os.path.join(dir_path, "%s_input.set.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的代谢集为文件，路径：%s" % (option_name, file_path))
    data = _get_objectid(data)
    main_id = db['metab_set'].find_one({"_id": data})["main_id"]
    main_id = _get_objectid(main_id)
    metab_set_list = db['metab_set_detail'].find_one({"set_id": main_id})['set_list']
    out_file = open(file_path,"w")
    out_file.write("metab_id" + "\n")
    for i in metab_set_list:
        out_file.write(i + "\n")
    out_file.close()
    return file_path


def export_paired_sample(data,option_name,dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, '%s_pair_sample.txt' % option_name)
    data_id = bind_obj.sheet.option('paired_id')
    collection = db['specimen_group_paired']
    search_r = collection.find_one({"_id": ObjectId(data_id)})
    if not search_r:
        raise Exception('specimen_group_paired中没有找到配对数据 %s'%data_id )
    with open(file_path ,'w') as fw:
        fw.write("#sample\tgroup\n")
        category_name = search_r['category_name']
        pair_category_name = search_r['pair_category_name']
        specimen_list = search_r['specimen_list'].split(',')
        pair_list = search_r['pair_list'].split(',')
        if len(specimen_list) != len(pair_list):
            raise Exception("specimen_list 和pair_list 长度不相等")
        for s in specimen_list:
            fw.write("%s\t%s\n"%(s,category_name))
        for s in pair_list:
            fw.write("%s\t%s\n"%(s,pair_category_name))
    return file_path


def export_asso_table(data, option_name, dir_path, bind_obj=None):

    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s表格为文件，路径:%s" % (option_name, file_path))

    group_detail = bind_obj.sheet.option('group_detail')
    if not isinstance(group_detail, dict):
        try:
            bind_obj.logger.debug("正确")
            table_dict = json.loads(group_detail)
            bind_obj.logger.debug(table_dict)
        except Exception:
            bind_obj.logger.debug("错误")
            raise Exception("传入的{}不是一个字典或者是字典对应的字符串".format('group_detail'))
    else:
        table_dict = group_detail

    assodata_id = ObjectId(data)
    assodata_main = db['assodata'].find_one({"_id": assodata_id})
    if not assodata_main:
        bind_obj.logger.debug("没有找到assodata %s"%data)
        raise Exception("没有找到assodata %s"%data)
    task_id = assodata_main['task_id']
    #var_list = assodata_main['assodata_names'].split(',')
    var_list = sorted(assodata_main['assodata_names'].keys(), key=lambda x: int(x.split('key')[1]))

    name_id_dic = dict()
    for s_name in table_dict.values():
        for s in s_name:
            ret_name = db['specimen'].find_one({'task_id':task_id, "name":s})
            if not ret_name:
                raise Exception('没有找到specimen 的样本名称%s'%(s))
            name_id_dic[s] = ret_name['_id']

    has_data_sample = []
    with open(file_path, "w") as f:
        f.write('sample\t'+'\t'.join(var_list)+'\n')
        for name in name_id_dic:
            tmp = [name]
            id = name_id_dic[name]
            assodata_detail = db['assodata_detail'].find_one({"assodata_id":assodata_id,"specimen_id":id})
            if not assodata_detail:
                #continue
                #bind_obj.set_error("关联表没有%s样本数据"%name)
                #bind_obj.logger.error("关联表没有%s样本数据"%name)
                raise Exception("关联表没有%s样本数据"%name)

            has_data_sample.append(name)
            for var in var_list:
                tmp.append(str(assodata_detail[var]))
            f.write('\t'.join(tmp)+'\n')
    if len(has_data_sample) <3:
        #bind_obj.set_error("有数据的样本数小于3")
        #bind_obj.logger.error("有数据的样本数小于3")
        raise Exception('有数据的样本数小于3')

    return file_path

#add by zhaoyuzhuo 2021.12.08
def export_hmdb_level(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s_input.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s为文件，路径：%s" % (option_name, file_path))
    f = open(file_path, 'w')
    main_id = db['anno_hmdb'].find_one({"task_id": data})["main_id"]
    result = db['anno_hmdb_level'].find({"hmdb_id": main_id})
    head_list = ['metab', 'metab_id', 'kingdom', 'superclass', 'class', 'subclass', 'o_id']
    f.write("\t".join(head_list) + "\n")
    for one in result:
        mongo_name_list = ['metab', 'metab_id', 'kingdom', 'superclass', 'class', 'subclass', 'o_id']
        write_list = map(lambda x: one[x], mongo_name_list)
        f.write("\t".join(write_list) + "\n")
    f.close()
    return file_path