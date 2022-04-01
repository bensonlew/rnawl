# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import re
import json
from types import StringTypes
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import defaultdict

client = Config().get_mongo_client(mtype="metaasv")
db = client[Config().get_mongo_dbname("metaasv")]
LEVEL = {
1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
6: "f__", 7: "g__", 8: "s__", 9: "asv"
}


def export_otu_table(data, option_name, dir_path, bind_obj=None):
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的OTU表格为文件，路径:%s" % (option_name, file_path))
    collection = db['asv_specimen']
    my_collection = db['specimen_detail']
    results = collection.find({"asv_id": ObjectId(data)})
    samples = []
    for result in results:
        id_result = my_collection.find_one({"_id": result["specimen_id"]})
        if not id_result:
            raise Exception("意外错误，样本id:{}在sg_specimen中未找到！")
        samples.append(id_result["specimen"])
    
    # samples = result["specimen_names"]
    # 因为有些样本名以1,2,3,4进行编号， 导致读出来了之后samples列表里的元素是数字， 需要先转化成字符串
    samples = map(str, samples)
    samples.sort()
    collection = db['asv_detail']
    with open(file_path, "wb") as f:
        f.write("ASV ID\t%s\ttaxonomy\n" % "\t".join(samples))
        for col in collection.find({"asv_id": ObjectId(data)}):
            line = "%s\t" % col["asv"]
            for s in samples:
                line += "%s\t" % col[s]
            for cls in ["d__", "k__", "p__", "c__", "o__", "f__", "g__"]:
                if cls in col.keys():
                    line += "%s; " % col[cls]
            line += col["s__"]
            f.write("%s\n" % line)
    client.close()
    return file_path


def export_otu_table_by_level(data, option_name, dir_path, bind_obj=None):
    """
    按等级获取OTU表
    使用时确保你的workflow的option里level这个字段
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的OTU表格为文件，路径:%s" % (option_name, file_path))
    collection = db['asv_specimen']
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_specimen表中未找到！".format(data))
    samples = list()
    for result in results:
        if "specimen_id" not in result:
            raise Exception("asv_id:{}错误，请使用新导入的ASV表的id".format(data))
        sp_id = result['specimen_id']
        my_collection = db['specimen_detail']
        my_result = my_collection.find_one({"_id": sp_id})
        if not my_result:
            raise Exception("意外错误，样本id:{}在specimen_detail表里未找到".format(sp_id))
        samples.append(my_result["specimen"])
        
    # 因为有些样本名以1,2,3,4进行编号， 导致读出来了之后samples列表里的元素是数字， 需要先转化成字符串
    samples = map(str, samples)
    samples.sort()
    level = int(bind_obj.sheet.option("level"))
    collection = db['asv_detail']
    name_dic = dict()
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_detail表中未找到！".format(data))
    for col in results:
        tmp = level + 1
        new_classify_name = _create_classify_name(col, tmp, bind_obj)
        if new_classify_name not in name_dic:
            name_dic[new_classify_name] = dict()
            for sp in samples:
                name_dic[new_classify_name][sp] = int(col[sp])
        else:
            for sp in samples:
                name_dic[new_classify_name][sp] += int(col[sp])
    with open(file_path, "wb") as f:
        f.write("ASV ID\t%s\n" % "\t".join(samples))
        for k in name_dic.iterkeys():
            line = k
            for s in samples:
                line += "\t" + str(name_dic[k][s])
            line += "\n"
            f.write(line)
    client.close()
    return file_path


def export_otu_table_by_level_remove_blank(data, option_name, dir_path, bind_obj=None):
    """
    OTU分析抽平用
    1.按等级获取OTU表，各层级以分号分割（不含有空格）
    ##fix by qingchen.zhang@20191120
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的ASV表格为文件，路径:%s" % (option_name, file_path))
    collection = db['asv_specimen']
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_specimen表中未找到！".format(data))
    samples = list()
    for result in results:
        if "specimen_id" not in result:
            raise Exception("asv_id:{}错误，请使用新导入的ASV表的id".format(data))
        sp_id = result['specimen_id']
        my_collection = db['specimen_detail']
        my_result = my_collection.find_one({"_id": sp_id})
        if not my_result:
            raise Exception("意外错误，样本id:{}在specimen_detail里未找到".format(sp_id))
        samples.append(my_result["specimen"])

    # 因为有些样本名以1,2,3,4进行编号， 导致读出来了之后samples列表里的元素是数字， 需要先转化成字符串
    samples = map(str, samples)
    samples.sort()
    level = int(bind_obj.sheet.option("level"))
    collection = db['asv_detail']
    choose_col = {s: 1 for s in samples}
    tax_col = []
    for k in range(1,10):
        if k <= level:
            choose_col[LEVEL[k]] = 1
            tax_col.append(LEVEL[k])
    results = collection.find({"asv_id": ObjectId(data)}, choose_col).batch_size(100000)
    with open(file_path, 'w') as w:
        w.write("ASV ID\t{}\n".format("\t".join(samples)))
        if level == 9:
            for col in results:
                asv_id = ";".join([col[tax] for tax in tax_col])
                w.write("{}\t{}\n".format(asv_id, "\t".join([str(col[sp]) for sp in samples])))
        else:
            asv_detail = {}
            for col in results:
                asv_id = ";".join([col[tax] for tax in tax_col])
                if asv_id not in asv_detail:
                    asv_detail[asv_id] = {sp: int(col[sp]) for sp in samples}
                else:
                    for sp in samples:
                        asv_detail[asv_id][sp] += int(col[sp])
            for asv_id in asv_detail:
                w.write("{}\t{}\n".format(asv_id, "\t".join([str(col[sp]) for sp in samples])))
    client.close()
    return file_path

def _create_classify_name_1(col, tmp, bind_obj):
    LEVEL = {
        1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
        6: "f__", 7: "g__", 8: "s__", 9: "asv"
    }
    for i in range(1, 10):
        if LEVEL[i] not in col:
            raise Exception("Mongo数据库中的taxonomy信息不完整")
    new_col = list()
    for i in range(1, tmp):
        new_col.append(col[LEVEL[i]])
    return ";".join(new_col) ##fix by qingchen.zhang@20191120

def _create_classify_name(col, tmp, bind_obj):
    LEVEL = {
        1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
        6: "f__", 7: "g__", 8: "s__", 9: "asv"
    }
    for i in range(1, 10):
        if LEVEL[i] not in col:
            raise Exception("Mongo数据库中的taxonomy信息不完整")
    new_col = list()
    for i in range(1, tmp):
        new_col.append(col[LEVEL[i]])
    return "; ".join(new_col)


def _get_only_classify_name(col, level, bind_obj):
    LEVEL = {
        1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
        6: "f__", 7: "g__", 8: "s__", 9: "asv"
    }
    if level in LEVEL:
        if LEVEL[level] in col:
            return col[LEVEL[level]]
        else:
            raise Exception('数据库中不存在列：{}'.format(LEVEL[level]))
    else:
        raise Exception('错误的分类水平：{}'.format(level))


def export_group_table(data, option_name, dir_path, bind_obj=None):
    """
    按group_id 和 组名获取group表
    使用时确保你的workflow的option里category_name这个字段
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    if data in ["all", "All", "ALL"]:
        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
        return file_path
    group_table = db['specimen_group']
    group_name_list = list()
    group_name = bind_obj.sheet.option("category_name")
    group_name_list = re.split(',', group_name)
    data = _get_objectid(data)
    group_schema = group_table.find_one({"_id": ObjectId(data)})
    if not group_schema:
        raise Exception("无法根据传入的group_id:{}在specimen_group表里找到相应的记录".format(data))
    c_name = group_schema["category_names"]

    with open(file_path, "wb") as f:
        f.write("#sample\t" + group_schema["group_name"] + "\n")

    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    index_list = _get_index_list(group_name_list, c_name)
    sample_id = list()  # [[id名,组名], [id名, 组名]，...]
    sample_name = list()  # [[样本名, 组名], [样本名, 组名], ...]
    for i in index_list:
        for k in group_schema['specimen_names'][i]:
            # k 样品ID
            sample_id.append([k, c_name[i]])
    for pair in sample_id:
        result = sample_table.find_one({"_id": ObjectId(pair[0])})
        if not result:
            raise Exception("无法根据传入的group_id:{}在样本表{}里找到相应的样本名".format(data, sample_table_name))
        sample_name.append([result['specimen'], pair[1]])
    sample_name.sort()
    my_data = bind_obj.sheet.data
    if "pan_id" in my_data['options']:
        sample_number_check(sample_name)

    with open(file_path, "ab") as f:
        for pair in sample_name:
            f.write("{}\t{}\n".format(pair[0], pair[1]))
    client.close()
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


def _get_index_list(group_name_list, c_name):
    """
    获取specimen_names字段的index (key)
    """
    length = len(c_name)
    index_list = list()
    for i in range(length):
        for g_name in group_name_list:
            if g_name == c_name[i]:
                index_list.append(i)
                break
    return index_list


def sample_number_check(sample_name):
    sample_count = defaultdict(int)
    for pair in sample_name:
        sample_count[pair[1]] += 1
    for k in sample_count:
        if sample_count[k] < 3:
            raise Exception("组{}里的样本数目小于三个，每个组里必须有三个以上的样本".format(k))


def export_group_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detal这个字段
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    if data in ["all", "All", "ALL"]:                #guanqing.zou 20180419
        group_detail = bind_obj.sheet.option('group_detail')
        table_dict = json.loads(group_detail)
        sample_table_name = 'specimen_detail'
        sample_table = db[sample_table_name]
        sample_name = list()

        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
            for k in table_dict:
                for sp_id in table_dict[k]:
                    sp = sample_table.find_one({"_id":ObjectId(sp_id)})
                    if not sp:
                        raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                    else:
                        sp_name = sp["specimen"]
                    sample_name.append(sp_name)
                    #f.write("{}\t{}\n".format(sp_name,sp_name))
            sample_name.sort()
            for specimen in sample_name:
                f.write("{}\t{}\n".format(specimen,"All"))
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

    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    sample_group = {}
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen"]
                sample_group[sp_name]=k
                #f.write("{}\t{}\n".format(sp_name, k))
        for key in sorted(sample_group.keys()):
            f.write("{}\t{}\n".format(key, sample_group[key]))
    client.close()
    return file_path

def export_group_table_by_detail_2(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detail这个字段和second_group_detail字段
    目的是生成一张group表(此处不做任何判断，完成信息的收集就可以了)
    该种格式的文件主要用于meta的sourcetracker分析
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')  #
    second_group_detail = bind_obj.sheet.option('second_group_detail')
    if not isinstance(group_detail, dict):
        try:
            table_dict_1 = json.loads(group_detail)
        except Exception:
            raise Exception("生成group1表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict_1, dict):
        raise Exception("生成group1表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(second_group_detail, dict):
        try:
            table_dict_2 = json.loads(second_group_detail)
        except Exception:
            raise Exception("生成group2表失败，传入的不是一个字典或者是字典对应的字符串")
    if not isinstance(table_dict_2, dict):
        raise Exception("生成group2表失败，传入的不是一个字典或者是字典对应的字符串")
    with open(file_path, "wb") as f:
        f.write("#SampleID\tEnv\tSourceSink\n")
    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict_1:
            for sp_id in table_dict_1[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen"]
                f.write("{}\t{}\tsource\n".format(sp_name, k))
    with open(file_path, "ab") as w:
        for k in table_dict_2:
            for sp_id_2 in table_dict_2[k]:
                sp_2 = sample_table.find_one({"_id": ObjectId(sp_id_2)})
                if not sp_2:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id_2, sample_table_name))
                else:
                    sp_name = sp_2["specimen"]
                w.write("{}\t{}\tsink\n".format(sp_name, k))
    return file_path

def export_group_table_by_detail_order(data, option_name, dir_path, bind_obj=None):
    """
    two_group分析wilcoxon符号秩和检验使用  add by qingchen.zhang @20190909
    1.根据参数group_detail样本顺序，从sg_specimen中获取样本名称，生成group_table表
    2.使用时确保workflow的option里有group_detail这个字段
    3.根据group_id检查sg_specimen_group有没有生成对应的表
    wilcoxon符号秩和检验要求样品配对使用，此方法生成的group表与页面分组方案中的样本顺序相同
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))

    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    group_table = db['specimen_group']
    table_dict = {}
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

    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen"]
                    f.write("{}\t{}\n".format(sp_name, k))
    client.close()
    return file_path

def export_cascading_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    根据group_detail生成group表或者二级group表
    使用时确保你的workflow的option里group_detail这个字段
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    second_group_detail = bind_obj.sheet.option('second_group_detail')
    group_table = db['specimen_group']
    if not isinstance(group_detail, dict):
        try:
            table_list = [json.loads(group_detail)]
        except Exception:
            raise Exception("生成group表失败，传入的一级分组不是一个字典或者是字典对应的字符串")
    if second_group_detail != '' and not isinstance(second_group_detail, dict):
        try:
            table_list.append(json.loads(second_group_detail))
        except Exception:
            raise Exception("生成group表失败，传入的二级分组不是一个字典或者是字典对应的字符串")
    # bind_obj.logger.debug('{}'.format(table_list))
    for i in table_list:
        if not isinstance(i, dict):
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    group_schema = group_table.find_one({"_id": ObjectId(data)})
    if not group_schema:
        raise Exception("无法根据传入的group_id:{}在specimen_group表里找到相应的记录".format(data))

    with open(file_path, "wb") as f:
        f.write("#sample\t" + group_schema["group_name"])

    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    _write_cascading_table(table_list, sample_table, file_path, sample_table_name)
    client.close()
    return file_path


def _write_class1_table(table_list, file_path, sample_table, sample_table_name):
    table_dict = table_list[0]
    with open(file_path, "ab") as f:
        f.write("\n")
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen"]
                f.write("{}\t{}\n".format(sp_name, k))


def _write_class2_table(table_list, file_path, sample_table, sample_table_name):
    dict1 = dict()  # 样本id -> 分组名
    dict2 = dict()
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
            sp = sample_table.find_one({"_id": ObjectId(k)})
            if not sp:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(k, sample_table_name))
            else:
                sp_name = sp["specimen"]
            f.write("{}\t{}\t{}\n".format(sp_name, dict1[k], dict2[k]))


def _write_cascading_table(table_list, sample_table, file_path, sample_table_name):
    length = len(table_list)
    if length == 1:
        _write_class1_table(table_list, file_path, sample_table, sample_table_name)
    elif length == 2:
        _write_class2_table(table_list, file_path, sample_table, sample_table_name)
    else:
        raise Exception("group_detal字段含有三个或以上的字典")


def export_otu_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按等级与分组信息(group_detail)获取OTU表
    使用时确保你的workflow的option里level与group_detail这个字段
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的OTU表格为文件，路径:%s" % (option_name, file_path))
    bind_obj.logger.debug(data)
    collection = db['asv_specimen']
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_specimen表中未找到！".format(data))
    samples = list()
    table_dict = {}
    group_detail = bind_obj.sheet.option("group_detail")
    bind_obj.logger.debug(group_detail)
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    sample_table = db['specimen_detail']
    for k in table_dict:
        for sp_id in table_dict[k]:
            sp = sample_table.find_one({"_id": ObjectId(sp_id)})
            if not sp:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, 'specimen_detail'))
            else:
                samples.append(sp["specimen"])
    samples.sort()

    level = int(bind_obj.sheet.option("level"))
    collection = db['asv_detail']
    name_dic = dict()
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_detail表中未找到！".format(data))
    for col in results:
        tmp = level + 1
        new_classify_name = _create_classify_name(col, tmp, bind_obj)
        if new_classify_name not in name_dic:
            name_dic[new_classify_name] = dict()
            for sp in samples:
                name_dic[new_classify_name][sp] = int(col[sp])
        else:
            for sp in samples:
                name_dic[new_classify_name][sp] += int(col[sp])
    with open(file_path, "wb") as f:
        f.write("ASV ID\t%s\n" % "\t".join(samples))
        for k in name_dic.iterkeys():
            line = k
            for s in samples:
                line += "\t" + str(name_dic[k][s])
            line += "\n"
            f.write(line)
    client.close()
    return file_path

def export_otu_table_by_detail2(data, option_name, dir_path, bind_obj=None):
    """
    按等级与分组信息(group_detail2)获取OTU表
    使用时确保你的workflow的option里level与group_detail2这个字段,为随机森林中“选择预测样品”参数生成丰度表（该丰度表中的样品与分组方案中样品必须没有交集）
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的OTU表格为文件，路径:%s" % (option_name, file_path))
    bind_obj.logger.debug(data)
    collection = db['asv_specimen']
    asv_id = bind_obj.sheet.option("asv_id")
    results = collection.find({"asv_id": ObjectId(asv_id)})
    if not results.count():
        raise Exception("asv_id: {}在asv_specimen表中未找到！".format(data))
    samples = list()
    table_dict = {}
    group_detail = bind_obj.sheet.option("predict_sample")
    bind_obj.logger.debug(group_detail)
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    sample_table = db['specimen_detail']
    for k in table_dict:
        for sp_id in table_dict[k]:
            sp = sample_table.find_one({"_id": ObjectId(sp_id)})
            if not sp:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, 'specimen_detail'))
            else:
                samples.append(sp["specimen"])
    samples.sort()

    level = int(bind_obj.sheet.option("level"))
    collection = db['asv_detail']
    name_dic = dict()
    results = collection.find({"asv_id": ObjectId(asv_id)})
    if not results.count():
        raise Exception("asv_id: {}在asv_detail表中未找到！".format(asv_id))
    for col in results:
        tmp = level + 1
        new_classify_name = _create_classify_name(col, tmp, bind_obj)
        if new_classify_name not in name_dic:
            name_dic[new_classify_name] = dict()
            for sp in samples:
                name_dic[new_classify_name][sp] = int(col[sp])
        else:
            for sp in samples:
                name_dic[new_classify_name][sp] += int(col[sp])
    with open(file_path, "wb") as f:
        f.write("ASV ID\t%s\n" % "\t".join(samples))
        for k in name_dic.iterkeys():
            line = k
            for s in samples:
                line += "\t" + str(name_dic[k][s])
            line += "\n"
            f.write(line)
    client.close()
    return file_path


def export_otu_table_by_detail_without_blank(data, option_name, dir_path, bind_obj=None):
    """
    按等级与分组信息(group_detail)获取OTU表
    使用时确保你的workflow的option里level与group_detail这个字段
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    ##去除最后的\t
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的OTU表格为文件，路径:%s" % (option_name, file_path))
    bind_obj.logger.debug(data)
    samples = list()
    table_dict = {}
    group_detail = bind_obj.sheet.option("group_detail")
    bind_obj.logger.debug(group_detail)
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    sample_table = db['specimen_detail']
    for k in table_dict:
        for sp_id in table_dict[k]:
            sp = sample_table.find_one({"_id": ObjectId(sp_id)})
            if not sp:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, 'specimen_detail'))
            else:
                samples.append(sp["specimen"])

    collection = db['asv_detail_level']
    name_dic = dict()
    results = collection.find({"asv_id": ObjectId(data), "level_id": 9})
    if not results.count():
        raise Exception("asv_id: {}在asv_detail表中未找到！".format(data))
    for col in results:
        tmp = 10
        new_classify_name = _create_classify_name(col, tmp, bind_obj)
        new_classify_name = new_classify_name.split(";")[-1].strip()
        if new_classify_name not in name_dic:
            name_dic[new_classify_name] = dict()
            name_dic[new_classify_name] = {sp : col[sp] for sp in samples}

    with open(file_path, "wb") as f:
        f.write("ASV ID\t%s\n" % "\t".join(samples))
        for col in name_dic.iterkeys():
            line_list = [col]+[str(name_dic[col][x]) for x in samples]
            new_line = "\t".join(line_list)
            f.write(new_line + "\n")
    client.close()
    return file_path


def export_otu_table_without_zero(data, option_name, dir_path, bind_obj=None):
    """
    按等级与分组信息(group_detail)获取OTU表
    使用时确保你的workflow的option里level与group_detail这个字段
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的OTU表格为文件，路径:%s" % (option_name, file_path))
    bind_obj.logger.debug(data)
    collection = db['asv_specimen']
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_specimen表中未找到！".format(data))
    samples = list()
    table_dict = {}
    group_detail = bind_obj.sheet.option("group_detail")
    bind_obj.logger.debug(group_detail)
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    sample_table = db['specimen_detail']
    for k in table_dict:
        for sp_id in table_dict[k]:
            sp = sample_table.find_one({"_id": ObjectId(sp_id)})
            if not sp:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, 'specimen_detail'))
            else:
                samples.append(sp["specimen"])
    samples.sort()

    level = int(bind_obj.sheet.option("level"))
    collection = db['asv_detail']
    name_dic = dict()
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_detail表中未找到！".format(data))
    for col in results:
        tmp = level + 1
        new_classify_name = _create_classify_name(col, tmp, bind_obj)
        if new_classify_name not in name_dic:
            name_dic[new_classify_name] = dict()
            for sp in samples:
                name_dic[new_classify_name][sp] = int(col[sp])
        else:
            for sp in samples:
                name_dic[new_classify_name][sp] += int(col[sp])
    with open(file_path, "wb") as f:
        f.write("ASV ID\t%s\n" % "\t".join(samples))
        for k in name_dic.iterkeys():
            line = k
            for s in samples:
                line += "\t" + str(name_dic[k][s])
            line += "\n"
            line_data = map(lambda x: float(x), line.strip("\n").split("\t")[1:])
            if not any(line_data):
                continue
            else:
                f.write(line)
    client.close()
    return file_path

def export_env_group_table(data, option_name, dir_path, bind_obj=None):
    """
    将分组信息与环境因子合并为一张表
    使用时确保你的workflow的option里group_detail这个字段,生成挑选样品的group文件
    确保存workflow的option里在group_detail_list或env_id与env_labs，用于生成环境因子表
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
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
        groups = []
        n = 0
        bind_obj.logger.info(group_detail_list)
        for group_id, group_detail in group_detail_list.items():
            # new_one_group = group_detail_list[one_group]
            one_group_dict = {}
            # for group_id, group_detail in new_one_group.items():
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
    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    with open(file_path, "wb") as f:
        f.write("#sample\t" + '\t'.join(head) + "\n")
        for sp_id in sample_list:
            sp = sample_table.find_one({"_id": ObjectId(sp_id)})
            if not sp:
                raise Exception("group_detail中的样本_id:{}在样本集中未找到".format(sp_id))
            else:
                sample_name = sp["specimen"]
            f.write(sample_name)
            try:
                bind_obj.sheet.option('group_detail_list')
                bind_obj.logger.info(groups)
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

def export_sample_list_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    根据group_detail生成样品list add by zouxuan 20180419
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
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
    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen"]
                if len(table_dict) == 1 and k in ['all', 'All', 'ALL']:
                    f.write("{}\t{}\n".format(sp_name, sp_name))
                else:
                    f.write("{}\t{}\n".format(sp_name, k))
    client.close()
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


def export_rep_seq_by_samples_abund(data, option_name,dir_path, bind_obj=None):
    LEVEL = {
        1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
        6: "f__", 7: "g__", 8: "s__", 9: "asv"
    }
    file_path = os.path.join(dir_path, "%s_rep.fasta" % option_name)

    group_detail = bind_obj.sheet.option('group_detail')
    table_dict = json.loads(group_detail)
    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    sample_name = list()

    for k in table_dict:
        for sp_id in table_dict[k]:
            sp = sample_table.find_one({"_id":ObjectId(sp_id)})
            if not sp:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
            else:
                sp_name = sp["specimen"]
                sample_name.append(sp_name)

    level = int(bind_obj.sheet.option("level"))
    level_key = LEVEL[level]

    collection = db['asv_detail']
    rep_seq = dict()

    otu_id = bind_obj.sheet.option('asv_id')
    for spe in set(data.split(',')):
        rep_seq[spe] = {'max_abund':0,'seq':''}
        results = collection.find({"asv_id": ObjectId(otu_id),level_key: spe})
        if not results.count():
            raise Exception("asv_id: {}在asv_detail表中未找到！".format(spe))
        for r in results:
            samples_abund = 0
            for sample in sample_name:
                samples_abund += float(r[sample])
            if samples_abund > rep_seq[spe]['max_abund']:
                rep_seq[spe]['max_abund'] = samples_abund
                rep_seq[spe]['seq'] = r['asv_rep']

    with open(file_path , 'w') as fw:
        for spe in rep_seq:
            fw.write('>%s\n%s\n'%(spe, rep_seq[spe]['seq']))

    return file_path


def get_seq_by_select_samples_pick_pre_nums(data, option_name, dir_path,bind_obj=None):
    LEVEL = {
        1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
        6: "f__", 7: "g__", 8: "s__", 9: "asv"
    }
    file_path = os.path.join(dir_path, "{}_rep.fasta".format(option_name))
    species_group =  os.path.join(dir_path, "species_group.xls")
    species_table = os.path.join(dir_path, "species_table.xls")
    group_detail = bind_obj.sheet.option('group_detail')
    table_dict = json.loads(group_detail)
    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    sample_name = dict()

    for k in table_dict:
        tmp_samples = []
        for sp_id in table_dict[k]:
            sp = sample_table.find_one({"_id":ObjectId(sp_id)})
            if not sp:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
            else:
                sp_name = sp["specimen"]
                tmp_samples.append(sp_name)
        sample_name[k] = tmp_samples

    otu_id = ObjectId(data)
    level_id = int(bind_obj.sheet.option("level"))
    level_k = LEVEL[level_id]
    #
    sg_otu_detail = db['asv_detail']
    ret_all = sg_otu_detail.find({"asv_id":otu_id})
    otu_sum = dict()
    sum_abund = dict()
    otu_seq = dict()
    ## 导出颜色组水平
    if int(bind_obj.sheet.option("color_level_id")) != 0:
        color_k = LEVEL[int(bind_obj.sheet.option("color_level_id"))]
        spe_color_group = dict()
    ## 物种在每组的丰度
    sum_abund_group = dict()
    for each in ret_all:
        spe = each[level_k]
        otu_name = each['asv']

        tmp_group_abund = {}
        for group in sample_name.keys():
            g_otu_abund = 0
            for sample in sample_name[group]:
                g_otu_abund += int(each[sample])
            tmp_group_abund[group] = g_otu_abund

        otu_abund = 0
        for g in tmp_group_abund:
            otu_abund += tmp_group_abund[g]

        if spe not in otu_sum:
            otu_sum[spe] = {otu_name:otu_abund}
            sum_abund[spe] = otu_abund
            if int(bind_obj.sheet.option("color_level_id")) != 0:
                spe_color_group[spe] = each[color_k]
            sum_abund_group[spe] = tmp_group_abund
        else:
            otu_sum[spe][otu_name] = otu_abund
            sum_abund[spe] += otu_abund
            for g in tmp_group_abund:
                sum_abund_group[spe][g] += tmp_group_abund[g]

        otu_seq[otu_name] = each['asv_rep']

    top_n = int(bind_obj.sheet.option("topN"))
    sum_abund_sort = sorted(sum_abund.iteritems(), key=lambda x: x[1],reverse=True)
    if len(sum_abund_sort) <= top_n: ##add by 2 line qingchen.zhang@20191209
        top_n = len(sum_abund_sort)
    top_n_list = []
    with open(file_path,'w') as fw:
        for i in range(top_n):
            pick_spe = sum_abund_sort[i][0]
            top_n_list.append(pick_spe)
            otu_abund_sort = sorted(otu_sum[pick_spe].iteritems(), key=lambda x : x[1], reverse=True)
            pick_otu_name = otu_abund_sort[0][0]
            seq = otu_seq[pick_otu_name]
            fw.write('>%s\n%s\n'%(pick_spe, seq))

    if int(bind_obj.sheet.option("color_level_id")) != 0:
        with open(species_group, 'w') as fw:
            fw.write("#species/ASV\tGROUP\n")
            for spe in top_n_list:
                fw.write(spe+"\t"+spe_color_group[spe]+"\n")

    with open(species_table, 'w') as fw:
        group_list = sample_name.keys()
        fw.write('ID\t'+'\t'.join(group_list)+"\n")
        for spe in top_n_list:
            tmp_line = [spe]
            for g in group_list:
                tmp_line.append(str(sum_abund_group[spe][g]))
            fw.write('\t'.join(tmp_line)+'\n')

    return file_path



def export_otu_seqs(data,  option_name,dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s_rep.fasta" % option_name)
    collection = db['asv_detail']
    with open(file_path, "wb") as f:
        for col in collection.find({"asv_id": ObjectId(data)}):
            f.write('>%s\n%s\n'%(col['asv'],col['asv_rep']))
    return file_path

def export_signal_pair_sample(data,option_name,dir_path, bind_obj=None):
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
    search_r = collection.find_one({"_id": ObjectId(data)})
    if not search_r:
        raise Exception('specimen_pair_group中没有找到配对数据 %s'%data )
    with open(file_path ,'w') as fw:
        category_name = search_r['category_name']
        pair_category_name = search_r['pair_category_name']
        fw.write('%s\t%s\n'%(category_name, pair_category_name))
        specimen_list = search_r['specimen_list'].split(',')
        pair_list = search_r['pair_list'].split(',')
        if len(specimen_list) != len(pair_list):
            raise Exception("specimen_list 和pair_list 长度不相等")
        spe_db= db['specimen_detail']
        for s, p in zip(specimen_list,pair_list):
            try:
                s_name = spe_db.find_one({"_id":ObjectId(s)})['specimen_name']
                p_name = spe_db.find_one({"_id":ObjectId(p)})['specimen_name']
            except Exception as e:
                raise  Exception("specimen_detail中没有找到%s或%s样本名称"%(s,p))
            fw.write("%s\t%s\n"%(s_name,p_name))
    return file_path

def export_est_table(data, option_name, dir_path, bind_obj=None):
    est_path = os.path.join(dir_path, "%s_input.estimators.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的多样性指数表格为文件，路径:%s" % (option_name, est_path))
    collection = db['alpha_diversity_detail']
    est_collection = db['alpha_diversity']
    task_id = "_".join(bind_obj.sheet.id.split('_')[0:2])
    result = est_collection.find_one({"_id": ObjectId(data), "task_id": task_id})
    if not result:
        raise Exception('没有找到多样性指数id对应的表，请检查传入的id是否正确')
    if not result['params']:
        index_type = u"ace,chao,shannon,simpson,coverage"
    elif type(result['params']) is dict:
        params = result["params"]
        if 'indices' in params:
            index_type = params['indices']
        elif 'index_type' in params:
            index_type = params['index_type']
    else:
        params = json.loads(result["params"])
        if 'indices' in params:
            index_type = params['indices']
        elif 'index_type' in params:
            index_type = params['index_type']
    indices = index_type.split(',')
    bind_obj.logger.debug(indices)
    details = collection.find({"alpha_diversity_id": ObjectId(data)})
    if not details.count():
        raise Exception('没有找到相应detail信息')
    index_nan = []
    write_data = {}
    for index in indices:
        if index == "jack":
            write_data["jackknife"] = ""
        else:
            write_data[index] = ""
    with open(est_path, "wb") as f:
        f.write("Estimators")
        for col in details:
            f.write("\t{}".format(col["specimen"]))
            for key in col:
                if key in write_data:
                    if str(col[key]) == "nan":
                        index_nan.append(key)
                    else:
                        write_data[key] += "\t{}".format(str(col[key]))
                else:
                    continue
        f.write("\n")
        for line in write_data:
            if line in index_nan:
                pass
            else:
                f.write("{}{}\n".format(line, write_data[line]))
    return est_path

def export_data_from_fastq_dir(data, option_name, dir_path, bind_obj=None):
    """
    根据前端传过来的fastq_dir转成file文件
    :param fastq_dir:
    :return:
    """
    fastq_dir_path = os.path.join(dir_path, "fastq.json")
    bind_obj.logger.debug("正在导出参数%s的file_list为文件夹，路径:%s" % (option_name, fastq_dir_path))
    with open(fastq_dir_path, 'w') as fw:
        fw.write(data)
    return fastq_dir_path

def export_bugbase_contribution_table(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的bugbase表格为文件，路径:%s" % (option_name, file_path))
    collection = db['bugbase']
    my_collection = db['bugbase_detail_contribution']
    results = my_collection.find({"bugbase_id": ObjectId(data)})

    with open(file_path, "wb") as f:
        f.write("#Query\tcategory\n")
        for col in results:
            phenotypes = col["phenotypes"]
            for xx in col:
                if col[xx] == 1:
                    f.write(str(xx) + "\t" + phenotypes + "\n")
    client.close()
    return file_path

def export_bugbase_table_by_bugbase_id(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的bugbase表格为文件，路径:%s" % (option_name, file_path))
    collection = db['bugbase']
    my_collection = db['bugbase_detail']
    samples = collection.find_one({"_id": ObjectId(data)})
    sample_list = str(samples["specimen_list"]).strip().split(",")
    del sample_list[0]
    results = my_collection.find({"bugbase_id": ObjectId(data)})

    with open(file_path, "wb") as f:
        f.write("Phenotypes\t" + "\t".join(sample_list) + "\n")
        for col in results:
            phenotypes = col["phenotypes"]
            f.write(phenotypes + "\t")
            tmp = []
            for xx in sample_list:
                tmp.append(str(col[xx]))
            f.write("\t".join(tmp) + "\n")
    client.close()
    return file_path

def export_faprotax_table_by_faprotax_id(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的bugbase表格为文件，路径:%s" % (option_name, file_path))
    collection = db['faprotax']
    my_collection = db['faprotax_detail']
    samples = collection.find_one({"_id": ObjectId(data)})
    sample_list = str(samples["specimen_list"]).strip().split(",")
    #del sample_list[0]
    results = my_collection.find({"faprotax_id": ObjectId(data)})

    with open(file_path, "wb") as f:
        f.write("Functional groups\t" + "\t".join(sample_list) + "\n")
        for col in results:
            if col["type"] == "sample":
                functional_groups = col["functional_groups"]
                f.write(functional_groups + "\t")
                tmp = []
                for xx in sample_list:
                    tmp.append(str(col[xx]))
                f.write("\t".join(tmp) + "\n")
    client.close()
    return file_path

def export_otu_table_without_tax(data, option_name, dir_path, bind_obj=None):
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的OTU表格为文件，路径:%s" % (option_name, file_path))
    collection = db['asv_specimen']
    my_collection = db['specimen_detail']
    results = collection.find({"asv_id": ObjectId(data)})
    samples = []
    for result in results:
        id_result = my_collection.find_one({"_id": result["specimen_id"]})
        if not id_result:
            raise Exception("意外错误，样本id:{}在sg_specimen中未找到！")
        samples.append(id_result["specimen"])

    # samples = result["specimen_names"]
    # 因为有些样本名以1,2,3,4进行编号， 导致读出来了之后samples列表里的元素是数字， 需要先转化成字符串
    samples = map(str, samples)
    samples.sort()
    collection = db['asv_detail']
    with open(file_path, "wb") as f:
        f.write("ASV ID\t%s\n" % "\t".join(samples))
        for col in collection.find({"asv_id": ObjectId(data)}):
            line = "%s\t" % col["asv"]
            for s in samples:
                line += "%s\t" % col[s]
            f.write("%s\n" % line)
    client.close()
    return file_path

def export_tax_table_by_asv_id(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的TAX表格为文件，路径:%s" % (option_name, file_path))
    collection = db['asv_detail']
    my_collection = db['specimen_detail']
    results = collection.find({"asv_id": ObjectId(data)})

    with open(file_path, "wb") as f:
        f.write("#Query\tDomain\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")
        for col in results:
            f.write(col["asv"] + "\t")
            for cls in ["d__", "k__", "p__", "c__", "o__", "f__", "g__"]:
                f.write(col[cls] +"\t")
            f.write(col["s__"] + "\n")
    client.close()
    return file_path

def export_asv_table_faprotax_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按等级与分组信息(group_detail)获取asv表
    使用时确保你的workflow的option里level与group_detail这个字段
    """
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的ASV表格为文件，路径:%s" % (option_name, file_path))
    bind_obj.logger.debug(data)
    collection = db['asv_specimen']
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_specimen表中未找到！".format(data))
    samples = list()
    table_dict = {}
    group_detail = bind_obj.sheet.option("group_detail")
    bind_obj.logger.debug(group_detail)
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    sample_table = db['specimen_detail']
    for k in table_dict:
        for sp_id in table_dict[k]:
            sp = sample_table.find_one({"_id": ObjectId(sp_id)})
            if not sp:
                raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, 'specimen_detail'))
            else:
                samples.append(sp["specimen"])
    samples.sort()

    level = 9
    collection = db['asv_detail']
    name_dic = dict()
    results = collection.find({"asv_id": ObjectId(data)})
    if not results.count():
        raise Exception("asv_id: {}在asv_detail表中未找到！".format(data))
    for col in results:
        tmp = level + 1
        new_classify_name = _create_classify_name(col, tmp, bind_obj)
        if new_classify_name not in name_dic:
            name_dic[new_classify_name] = dict()
            for sp in samples:
                name_dic[new_classify_name][sp] = int(col[sp])
        else:
            for sp in samples:
                name_dic[new_classify_name][sp] += int(col[sp])
    with open(file_path, "wb") as f:
        f.write("#ASV ID\t%s\ttaxonomy\n" % "\t".join(samples))
        for k in name_dic.iterkeys():
            line = k.split(";")[-1].strip(" ")
            for s in samples:
                line += "\t" + str(name_dic[k][s])
            line += "\t" + ";".join(k.strip().split(";")[0:-1])
            line += "\n"
            f.write(line)
    client.close()
    return file_path

def export_group_table_by_detail_unsort(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detal这个字段
    """
    # client = Config().mongo_client
    # db = client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    if data in ["all", "All", "ALL"]:                #guanqing.zou 20180419
        group_detail = bind_obj.sheet.option('group_detail')
        table_dict = json.loads(group_detail)
        sample_table_name = 'specimen_detail'
        sample_table = db[sample_table_name]
        sample_name = list()

        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
            for k in table_dict:
                for sp_id in table_dict[k]:
                    sp = sample_table.find_one({"_id":ObjectId(sp_id)})
                    if not sp:
                        raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                    else:
                        sp_name = sp["specimen"]
                    sample_name.append(sp_name)
                    #f.write("{}\t{}\n".format(sp_name,sp_name))
            sample_name.sort()
            for specimen in sample_name:
                f.write("{}\t{}\n".format(specimen,"All"))
        return file_path

    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    sample_list = []
    for i in group_detail.split(":")[:-1]:
        if i.strip():
            sample = i.strip().replace("{","").split(",")[-1].replace('"','')
            sample_list.append(sample)
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

    sample_table_name = 'specimen_detail'
    sample_table = db[sample_table_name]
    sample_group = {}
    with open(file_path, "ab") as f:
        for k in sample_list:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen"]
                sample_group[sp_name]=k
                f.write("{}\t{}\n".format(sp_name, k))
        #for key in sorted(sample_group.keys()):
        #    f.write("{}\t{}\n".format(key, sample_group[key]))
    client.close()
    return file_path

def export_normalized_table_by_bugbase_id(data, option_name, dir_path, bind_obj=None):
    """
    按bugbase_id获取bugbase normalized表
    使用时确保你的workflow的option里bugbase_id这个字段
    """
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的normalized表格为文件，路径:%s" % (option_name, file_path))
    collection = db['bugbase']
    main_rsult = collection.find_one({"main_id": ObjectId(data)})
    my_collection = db['bugbase_detail_normalized']
    results = my_collection.find({"bugbase_id": ObjectId(data)})
    samples = []
    specimen_list = str(main_rsult["specimen_list"])
    for i in specimen_list.split(","):
        if i != "phenotypes":
            samples.append(i)

    with open(file_path, "wb") as f:
        f.write("otu id" + "\t" + "\t".join(samples) + "\n")
        for col in results:
            f.write(col["otu_id"])
            for xx in samples:
                f.write("\t" + col[xx])
            f.write("\n")
    client.close()
    return file_path