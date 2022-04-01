# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last modified BY LIULinmeng, for mongodb, 20171113
import os
import re
import json
from types import StringTypes
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import defaultdict
from mainapp.models.mongo.metagenomic import Metagenomic
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mbio.packages.metagenomic.id_convert import id2name, get_group_collection

client = Config().get_mongo_client(mtype="metagenomic")
db = client[Config().get_mongo_dbname("metagenomic")]


def export_group_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detail这个字段
    data 为task_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件 add by zhujuan 2017.11.09
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB + '_metagenomic']
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    metagenomic = Metagenomic()
    metagenomic._config = Config()
    geneset_info = metagenomic.get_geneset_info(data)
    task_id = geneset_info['task_id']
    id_to_sample = id2name(task_id, type="task")
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
    with open(file_path, "ab") as f:
        for k in sorted(table_dict):
            for sp_id in table_dict[k]:
                sample_name = id_to_sample[sp_id]  # id2name函数输出的id类型变换了
                if not sample_name:
                    raise Exception("group_detail中的样本_id:{}在样本集中未找到".format(sp_id))
                if len(table_dict) == 1 and k in ['all', 'All', 'ALL']:
                    f.write("{}\t{}\n".format(sample_name, sample_name))
                else:
                    f.write("{}\t{}\n".format(sample_name, k))
    client.close()
    return file_path
    """
    if data in ["all", "All", "ALL"]:
        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
        return file_path
    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    group_table = db['specimen_group']
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    group_schema = group_table.find_one({"_id": ObjectId(data)})
    if not group_schema:
        raise Exception("无法根据传入的group_id:{}在specimen_group表里找到相应的记录".format(data))
    schema_name = re.sub("\s", "_", group_schema["group_name"])  # 将分组方案名的空格替换成下划线
    with open(file_path, "wb") as f:
        f.write("#sample\t" + schema_name + "\n")

    sample_table_name = 'data_stat_detail'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    raise Exception("group_detail中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen_name"]
                f.write("{}\t{}\n".format(sp_name, k))
    client.close()
    return file_path
    """


def export_group_table_for_signal(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表,符号秩和检验(signal)为配对检验，对样本顺序有要求，所以按group表中的样本顺序进行排序
    使用时确保你的workflow的option里group_detail这个字段
    data 为task_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件 add by guhaidong 2017.12.15
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB + '_metagenomic']
    geneset_id, group_id = data.split(":")
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在运行按分组样品顺序的导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data = _get_objectid(geneset_id)
    group_id = _get_objectid(group_id)
    group_detail = bind_obj.sheet.option('group_detail')
    metagenomic = Metagenomic()
    metagenomic._config = Config()
    geneset_info = metagenomic.get_geneset_info(data)
    task_id = geneset_info['task_id']
    id_to_sample = id2name(task_id, type="task")
    tmp, sp_list_tmp = get_group_collection(group_id)
    sp_list = []
    for sub_list in sp_list_tmp:
        sp_list += sub_list
    # bind_obj.logger.debug("sp_list debug print:")
    # bind_obj.logger.debug(sp_list)
    if not isinstance(group_detail, dict):
        try:
            bind_obj.logger.debug("group_detail格式正确")
            table_dict = json.loads(group_detail)
        except Exception:
            bind_obj.logger.debug("group_detail格式错误")
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    else:
        table_dict = group_detail
    if not isinstance(table_dict, dict):
        bind_obj.logger.debug("table_dict格式错误")
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    schema_name = "group_name"
    with open(file_path, "wb") as f:
        f.write("#sample\t" + schema_name + "\n")
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in sp_list:
                if sp_id in table_dict[k]:
                    sample_name = id_to_sample[sp_id]
                    f.write("{}\t{}\n".format(sample_name, k))

    client.close()
    return file_path


def export_env_group_table(data, option_name, dir_path, bind_obj=None):
    """
    将分组信息与环境因子合并为一张表
    使用时确保你的workflow的option里group_detail这个字段,生成挑选样品的group文件
    确保存workflow的option里在group_detail_list或env_id与env_labs，用于生成环境因子表
    data 为task_id ;
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB + '_metagenomic']
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    data = _get_objectid(data)
    metagenomic = Metagenomic()
    metagenomic._config = Config()
    geneset_info = metagenomic.get_geneset_info(data)
    task_id = geneset_info['task_id']
    id_to_sample = id2name(task_id, type="task")
    head = []
    sample_dict = bind_obj.sheet.option('group_detail')
    if not isinstance(sample_dict, dict):
        try:
            bind_obj.logger.debug("正确")
            sample_dict = json.loads(sample_dict)
            bind_obj.logger.debug(sample_dict)
        except Exception:
            bind_obj.logger.debug("错误")
            raise Exception("传入的group_detail不是一个字典或者是字典对应的字符串")
    if not isinstance(sample_dict, dict):
        bind_obj.logger.debug("错误")
        raise Exception("生成的group_detail不是一个字典或者是字典对应的字符串")
    if sample_dict.keys()[0] not in ['all', 'All', 'ALL']:
        raise Exception("只能在all下选择的样品")
    else:
        sample_list = sample_dict.values()[0]
    try:
        bind_obj.sheet.option('group_detail_list')
        bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
        group_detail_list = eval(bind_obj.sheet.option('group_detail_list'))
        # group_list = bind_obj.sheet.option('group_list')
        group_collection = db['specimen_group']
        group_dict = {}
        # for group_id in group_list.split(","):
        #     group_info = group_collection.find_one({'_id': _get_objectid(group_id)})
        #     group_name.append(str(group_info['group_name']))
        groups=[]
        n = 0
        bind_obj.logger.info(group_detail_list)
        for one_group in group_detail_list:
            one_group_dict = {}
            for group_id, group_detail in one_group.items():
                group_info = group_collection.find_one({'_id': _get_objectid(group_id)})
                group_name = str(group_info['group_name'])
                for i in group_detail.keys():
                    for j in group_detail[i]:
                        if group_name in ['all', 'All', 'ALL']:
                            one_group_dict[j] = j
                        else:
                            one_group_dict[j] = i
                group_dict[group_name] = one_group_dict
                head.append(group_name)
                groups.append(group_name)
    except:
        pass
    try:
        bind_obj.sheet.option('env_id')
        env_id = bind_obj.sheet.option('env_id')
        env_labs = bind_obj.sheet.option("env_labs").split(',')
        env_collection = db['env_detail']
        head.extend(env_labs)
    except:
        pass
    with open(file_path, "wb") as f:
        f.write("#sample\t" + '\t'.join(head) + "\n")
        for sp_id in sample_list:
            sample_name = id_to_sample[sp_id]  # id2name函数输出的id类型变换了
            # bind_obj.logger.info(sample_name)
            if not sample_name:
                raise Exception("group_detail中的样本_id:{}在样本集中未找到".format(sp_id))
            f.write(sample_name)
            try:
                bind_obj.sheet.option('group_detail_list')
                for k in groups:
                    f.write("\t" + group_dict[k][sp_id])
            except:
                bind_obj.logger.info("不存在分组相关信息")
            try:
                bind_obj.sheet.option('env_id')
                for e in env_labs:
                    result = env_collection.find_one({'env_id': _get_objectid(env_id),
                                                      'specimen_id': _get_objectid(sp_id)})
                    f.write("\t" + str(result[e]))
            except:
                bind_obj.logger.info("不存在环境因子相关信息")
            f.write("\n")
    client.close()
    return file_path


def export_env_table(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + '_metagenomic']
    file_path = os.path.join(dir_path, '%s_input_env.xls' % option_name)
    bind_obj.logger.debug('正在导出参数%s的环境因子表为文件，路径:%s' % (option_name, file_path))
    collection_main = db['env']
    result_main = collection_main.find_one({'_id': _get_objectid(data)})
    if not result_main:
        raise Exception('环境因子id没有找到对应的表信息')
    all_envs = result_main['env_names'].strip().split(',')
    collection = db['env_detail']
    specimen_collection = db['data_stat_detail']
    results = collection.find({'env_id': _get_objectid(data)})
    if results.count() == 0:
        raise Exception('环境因子id没有找到对应的detail数据')
    bind_obj.logger.info('ALL ENVS:' + ' '.join(all_envs))
    with open(file_path, 'wb') as f:
        f.write('#SampleID\t' + '\t'.join(all_envs) + '\n')
        for one in results:
            specimen_name = specimen_collection.find_one({'_id': _get_objectid(one['specimen_id'])})
            if not specimen_name:
                specimen_name = specimen_collection.find_one({'origin_id': _get_objectid(one['specimen_id'])})
            if specimen_name:
                line_list = [specimen_name['specimen_name']]
                for env in all_envs:
                    line_list.append(str(one[env]))
                f.write('\t'.join(line_list) + '\n')
            else:
                raise Exception('样本id:%s在data_stat_detail表中没有找到' % str(one['specimen_id']), )
    return file_path


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


def export_float_env(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + '_metagenomic']
    file_path = os.path.join(dir_path, '%s_input_env.xls' % option_name)
    bind_obj.logger.debug('正在导出参数%s的环境因子表为文件，路径:%s' % (option_name, file_path))
    collection_main = db['env']

    bind_obj.logger.debug(bind_obj.sheet.option("env_labs"))
    env_labs = bind_obj.sheet.option("env_labs").split(',')

    result_main = collection_main.find_one({'_id': _get_objectid(data)})
    if not result_main:
        raise Exception('环境因子id没有找到对应的表信息')
    all_envs = result_main['env_names'].strip().split(',')
    collection = db['env_detail']
    specimen_collection = db['data_stat_detail']
    results = collection.find({'env_id': _get_objectid(data)})
    if results.count() == 0:
        raise Exception('环境因子id没有找到对应的detail数据')
    bind_obj.logger.info('ALL ENVS:' + ' '.join(all_envs))
    bind_obj.logger.info('SELECT ENVS:' + ' '.join(env_labs))
    # write_lines = ['#SampleID\t' + '\t'.join(env_labs) + '\n']
    # flit_envs = set()
    with open(file_path, 'wb') as f:
        f.write('#SampleID\t' + '\t'.join(env_labs) + '\n')
        for one in results:
            specimen_name = specimen_collection.find_one({'_id': _get_objectid(one['specimen_id'])})
            if not specimen_name:
                specimen_name = specimen_collection.find_one({'origin_id': _get_objectid(one['specimen_id'])})
            if specimen_name:
                line_list = [specimen_name['specimen_name']]
                for env in env_labs:
                    line_list.append(str(one[env]))
                f.write('\t'.join(line_list) + '\n')
            else:
                raise Exception('样本id:%s在data_stat_detail表中没有找到' % str(one['specimen_id']))
    return file_path


def export_abund_table_path(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + '_metagenomic']
    bind_obj.logger.debug('正在导出参数%s的丰度文件路径' % option_name)
    collection_main = db['abund_table_path']
    result_main = collection_main.find_one({'_id': _get_objectid(data)})
    if not result_main:
        raise Exception('丰度文件没有找到对应的表信息')
    abund_table_path = result_main['abu_file']
    if not abund_table_path:
        raise Exception('找不到abund_table文件')
    return abund_table_path


def export_cascading_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    根据group_detail生成group表或者二级group表
    使用时确保你的workflow的option里group_detail这个字段
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data2 = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    group_id = bind_obj.sheet.option('group_id')
    second_group_detail = bind_obj.sheet.option('second_group_detail')
    group_table = db['specimen_group']
    table_list = ''
    if not isinstance(group_detail, dict):
        try:
            table_list = [json.loads(group_detail)]
        except Exception:
            raise Exception("生成group表失败，传入的一级分组不是一个字典或者是字典对应的字符串")
    if second_group_detail != '"null"' and not isinstance(second_group_detail, dict):
        try:
            table_list.append(json.loads(second_group_detail))
        except Exception:
            raise Exception("生成group表失败，传入的二级分组不是一个字典或者是字典对应的字符串")
    # bind_obj.logger.debug('{}'.format(table_list))
    for i in table_list:
        if not isinstance(i, dict):
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    group_schema = group_table.find_one({"_id": _get_objectid(group_id)})
    if not group_schema:
        raise Exception("无法根据传入的group_id:{}在specimen_group表里找到相应的记录".format(data2))

    with open(file_path, "wb") as f:
        f.write("#sample\t" + group_schema["group_name"])
    # bind_obj.logger.debug('{}'.format(table_list))
    _write_cascading_table(data2, table_list, file_path)
    client.close()
    return file_path


def _write_class1_table(data, table_list, file_path):
    table_dict = table_list[0]
    data2 = _get_objectid(data)
    metagenomic = Metagenomic()
    metagenomic._config = Config()
    geneset_info = metagenomic.get_geneset_info(data2)
    task_id = geneset_info['task_id']
    id_to_sample = id2name(task_id, type="task")
    with open(file_path, "ab") as f:
        f.write("\n")
        for k in table_dict:
            for sp_id in table_dict[k]:
                sample_name = id_to_sample[sp_id]  # id2name函数输出的id类型变换了
                if not sample_name:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sample_name, k))
                else:
                    sample_name = sample_name
                f.write("{}\t{}\n".format(sample_name, k))


def _write_class2_table(data, table_list, file_path):
    dict1 = dict()  # 样本id -> 分组名
    dict2 = dict()
    data2 = _get_objectid(data)
    metagenomic = Metagenomic()
    metagenomic._config = Config()
    geneset_info = metagenomic.get_geneset_info(data2)
    task_id = geneset_info['task_id']
    id_to_sample = id2name(task_id, type="task")
    for k in table_list[0]:
        for id_ in table_list[0][k]:
            if id_ in dict1:
                raise Exception("样本id:{}在分组{}，{}中都出现，一个样本在同一级别上只能属于一个分组！".format(id_, dict1[id_], k))
            dict1[id_] = k
    for k in table_list[1]:
        for id_ in table_list[1][k]:
            if id_ not in dict1:
                raise Exception("样本id:{}在第一级的分组中未出现".format(id_))
            if id_ in dict2:
                raise Exception("样本id:{}在分组{}，{}中都出现，一个样本在同一级别上只能属于一个分组！".format(id_, dict1[id_], k))
            dict2[id_] = k
    if len(dict1) != len(dict2):
        raise Exception("一级分组中的某些样本id在二级分组中未找到！")
    with open(file_path, "ab") as f:
        f.write("\tsecond_group\n")
        for k in dict1:
            sample_name = id_to_sample[k]
            if not sample_name:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sample_name, k))
            else:
                sample_name = sample_name
            f.write("{}\t{}\t{}\n".format(sample_name, dict1[k], dict2[k]))


def _write_cascading_table(data, table_list, file_path):
    length = len(table_list)
    if length == 1:
        _write_class1_table(data, table_list, file_path)
    elif length == 2:
        _write_class2_table(data, table_list, file_path)
    else:
        raise Exception("group_detal字段含有三个或以上的字典")

def export_function_set(data, option_name, dir_path, bind_obj=None):
    """
    获得个性化功能集表
    data: {"function_id": 5be4f0b3d4daa48652b74365, "type":"1"}
    """
    #data = json.loads(data)
    data = json.loads(data)[0]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的代谢集为文件，路径：%s" % (option_name, file_path))
    function_id = data["function_id"]
    level_type = data["type"]
    function_id = _get_objectid(function_id)
    result = db['func_set'].find_one({"_id": function_id})
    if not result:
        raise Exception("func_set not find main table : {}".format(function_id))
    print function_id
    anno_type = result["anno_type"]
    if anno_type == "cog":
        lowest_level = "nog_list"
        level = "NOG"
    elif anno_type == "kegg":
        lowest_level = "k_list"
        #level = "Gene"
        level = "KO"
    elif anno_type == "cazy":
        lowest_level = "family_list"
        level = "Family"
    elif anno_type == "ardb":
        lowest_level = "arg_list"
        level = "ARG"
    elif anno_type == "card":
        lowest_level = "aro_list"
        level = "ARO"
    elif anno_type == "vfdb":
        lowest_level = "vfs_list"
        level = "VFs"
    if not result[lowest_level]:
        raise Exception("func_set_id -{} has no {} arg".format(function_id, lowest_level))
    func_set_all = result[lowest_level]
    out_file = open(file_path,"w")
    for name in func_set_all:
        out_file.write(level_type + "\t" + level + "\t" + name + "\n")
    out_file.close()
    return file_path

def export_group_table_by_rf_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detail这个字段
    data 为task_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件 add by zhujuan 2017.11.09
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB + '_metagenomic']
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('pre_group_detail')
    metagenomic = Metagenomic()
    metagenomic._config = Config()
    geneset_info = metagenomic.get_geneset_info(data)
    task_id = geneset_info['task_id']
    id_to_sample = id2name(task_id, type="task")
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
    with open(file_path, "ab") as f:
        for k in sorted(table_dict):
            for sp_id in table_dict[k]:
                sample_name = id_to_sample[sp_id]  # id2name函数输出的id类型变换了
                if not sample_name:
                    raise Exception("group_detail中的样本_id:{}在样本集中未找到".format(sp_id))
                if len(table_dict) == 1 and k in ['all', 'All', 'ALL']:
                    f.write("{}\t{}\n".format(sample_name, sample_name))
                else:
                    f.write("{}\t{}\n".format(sample_name, k))
    client.close()
    return file_path

def export_clean_data(data, option_name, dir_path, bind_obj=None):
    """
    获得clean_stat
    data: task_id
    """
    task_id = data
    file_path = os.path.join(dir_path, "%s__%s.xls" % (option_name, task_id))
    bind_obj.logger.debug("正在导出参数%s的代谢集为文件，路径：%s" % (option_name, file_path))
    result = db['data_stat'].find_one({"task_id": task_id, "type" : "clean"})
    if not result:
        raise Exception("func_set not find main table : {}".format(function_id))
    detail_id = result["_id"]
    detail = db['data_stat_detail'].find({"data_stat_id": detail_id})
    if not detail:
        raise Exception("not find main detail : {}".format(detail_id))
    samples_dic = id2name(task_id, type="task")
    out_file = open(file_path,"w")
    out_file.write("#Sample\tReadsNum\n")
    for each in detail:
        specimen_id = each["specimen_name"]
        clean_read_num = each["clean_read_num"]
        if samples_dic.has_key(specimen_id):
            specimen_name = samples_dic[specimen_id]
        else:
            raise Exception("specimen_id : {} not in samples_dic".format(specimen_id))
        out_file.write(specimen_name + "\t" + str(clean_read_num) + "\n")
    out_file.close()
    return file_path


def export_env_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取env 分组表
    使用时确保你的workflow的option里env_detail这个字段
    zouguanqing
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的env group表格为文件，路径:%s" % (option_name, file_path))

    group_detail = bind_obj.sheet.option('env_detail')
    #metagenomic = Metagenomic()
    #geneset_info = metagenomic.get_geneset_info(data)
    #task_id = geneset_info['task_id']
    #id_to_sample = id2name(task_id, type="task")
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

    with open(file_path, "wb") as f:
        f.write("#env\tgroup\n")
        for k in sorted(table_dict):
            for env in table_dict[k]:
                    f.write("{}\t{}\n".format(env, k))
    return file_path


def export_env_table_select(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + '_metagenomic']
    # workflow must has self.option('env_labs')
    file_path = os.path.join(dir_path, '%s_input_env.xls' % option_name)
    bind_obj.logger.debug('正在导出参数%s的环境因子表为文件，路径:%s' % (option_name, file_path))
    collection_main = db['env']
    result_main = collection_main.find_one({'_id': _get_objectid(data)})
    if not result_main:
        raise Exception('环境因子id没有找到对应的表信息')
    #all_envs = result_main['env_names'].strip().split(',')
    all_envs = bind_obj.sheet.option('env_labs').split(',')
    collection = db['env_detail']
    specimen_collection = db['data_stat_detail']
    results = collection.find({'env_id': _get_objectid(data)})
    if results.count() == 0:
        raise Exception('环境因子id没有找到对应的detail数据')
    bind_obj.logger.info('Select ENVS:' + ' '.join(all_envs))
    with open(file_path, 'wb') as f:
        f.write('#SampleID\t' + '\t'.join(all_envs) + '\n')
        for one in results:
            specimen_name = specimen_collection.find_one({'_id': _get_objectid(one['specimen_id'])})
            if not specimen_name:
                specimen_name = specimen_collection.find_one({'origin_id': _get_objectid(one['specimen_id'])})
            if specimen_name:
                line_list = [specimen_name['specimen_name']]
                for env in all_envs:
                    line_list.append(str(one[env]))
                f.write('\t'.join(line_list) + '\n')
            else:
                raise Exception('样本id:%s在data_stat_detail表中没有找到' % str(one['specimen_id']), )
    return file_path

def export_group_table_by_detail_2(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detail这个字段
    data 为task_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件 add by zhujuan 2017.11.09
    all时 分为1组
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB + '_metagenomic']
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    metagenomic = Metagenomic()
    metagenomic._config = Config()
    geneset_info = metagenomic.get_geneset_info(data)
    task_id = geneset_info['task_id']
    id_to_sample = id2name(task_id, type="task")
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
    with open(file_path, "ab") as f:
        for k in sorted(table_dict):
            for sp_id in table_dict[k]:
                sample_name = id_to_sample[sp_id]  # id2name函数输出的id类型变换了
                if not sample_name:
                    raise Exception("group_detail中的样本_id:{}在样本集中未找到".format(sp_id))
                if len(table_dict) == 1 and k in ['all', 'All', 'ALL']:
                    f.write("{}\t{}\n".format(sample_name, k))   #  此处有别于export_group_table_by_detail
                else:
                    f.write("{}\t{}\n".format(sample_name, k))
    client.close()
    return file_path


def get_mapping_file(data, option_name, dir_path, bind_obj=None):
    import urllib
    mapping_file = os.path.join(dir_path, "mapping_file.txt")
    mp_str = "sequence.fastq_dir||filelist[in_fastq](sanger):{}{{Rawdata}}".format(mapping_file)
    bind_obj.logger.info("正在设置%s的mapping_file路径:%s" % (option_name, mp_str))
    signature = {
        "ceshi": 1,
        "task_id": bind_obj.sheet.option('task_id'),
        "params_path": data
    }
    signature = urllib.urlencode(signature)
    if "sanger-dev" in os.path.abspath(__file__):
        url = "http://openapi.nsg.com/file/check_mapping_file?{}".format(signature)
    else:
        url = "http://openapi.labsanger.sanger.com/file/check_mapping_file?{}".format(signature)
    try:
        resp = urllib.urlopen(url)
    except Exception as e:
        raise Exception('下载mapping_file失败,{}'.format(e))

    contents = json.load(resp)
    if contents["success"] not in ['true', True]:
        raise Exception('下载mapping_file失败')
    with open(mapping_file, 'w') as w:
        json_idc = {option_name: contents['d']}
        json.dump(json_idc, w)
    return mp_str


def export_group_by_detail_taxon(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    name_id = bind_obj.sheet.option("name2id")
    id_name = {v: k for k, v in json.loads(name_id).items()}
    group_detail = json.loads(bind_obj.sheet.option("group"))
    with open(file_path, 'w') as w:
        w.write("#sample\tgroup_name\n")
        for k, v in group_detail.items():
            for i in v:
                if k.lower() == "all":
                    k = i
                w.write("{}\t{}\n".format(id_name[i], k))
    return file_path


def export_signal_pair_sample(data, option_name, dir_path, bind_obj=None):
    """
    配对样本的记录
    :param data:
    :param option_name:
    :param dir_path:
    :param bind_obj:
    :return:
    """
    file_path = os.path.join(dir_path, '%s_pair_sample.txt' % option_name)
    collection = db['specimen_pair_group']
    name2id = json.loads(bind_obj.sheet.option("name2id"))
    id2name = {i: n for n, i in name2id.items()}
    search_r = collection.find_one({"_id": ObjectId(data)})
    if not search_r:
        raise Exception('specimen_pair_group中没有找到配对数据 %s'%data )
    with open(file_path, 'w') as fw:
        category_name = search_r['category_name']
        pair_category_name = search_r['pair_category_name']
        fw.write('%s\t%s\n' % (category_name, pair_category_name))
        specimen_list = search_r['specimen_list'].split(',')
        pair_list = search_r['pair_list'].split(',')
        if len(specimen_list) != len(pair_list):
            raise Exception("specimen_list 和pair_list 长度不相等")
        for s, p in zip(specimen_list, pair_list):
            try:
                s_name = id2name[s]
                p_name = id2name[p]
            except Exception as e:
                raise Exception("specimen_detail中没有找到%s或%s样本名称"%(s, p))
            fw.write("%s\t%s\n" % (s_name, p_name))
    return file_path
