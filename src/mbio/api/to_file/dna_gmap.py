# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last modify 20180705

import os
import json
from bson.objectid import ObjectId
from biocluster.config import Config


db = Config().get_mongo_client(mtype="dna_gmap")[Config().get_mongo_dbname("dna_gmap")]


def export_trit_file(data, option_name, dir_path, bind_object=None):
    trit_coll = db["sg_feature_file"]
    trit_detail_coll = db["sg_feature_file_detail"]
    trit_file = os.path.join(dir_path, "trait.txt")
    result = trit_coll.find_one({"_id": ObjectId(data)})
    if not result:
        raise Exception("没有在表sg_feature_file中找到id:{}对应的结果".format(data))
    feature_dict = result["feature_dict"]
    specimen_ids = result["specimen_list"]
    task_id = result["task_id"]
    trit_list, value_list = [], []
    for trit in feature_dict.keys():
        trit_list.append(trit)
        value_list.append(feature_dict[trit])
    specimen_list = []
    child_coll = db["sg_specimen_child"]
    results = child_coll.find({"task_id": task_id})
    if not results:
        raise Exception("没有在表sg_specimen_child中找到task_id:{}对应的结果".format(task_id))
    for result in results:
        analysis_name = result["analysis_name"]
        if analysis_name in specimen_ids:
            specimen_list.append(analysis_name)
    specimen_list = list(set(specimen_list))
    with open(trit_file, "w") as w:
        w.write("#SampleID\t" + "\t".join(value_list) + "\n")
        # for sample in specimen_ids:
        for sample in specimen_list:
            result = trit_detail_coll.find_one({"feature_file_id": ObjectId(data), "sample": sample})
            if not result:
                raise Exception("没有在表sg_feature_file_detail表中找到样本：{}对应的结果".format(sample))
            line = []
            line.append(result["sample"])
            for t in trit_list:
                try:
                    line.append(result[t])
                except:
                    line.append("NaN")
            w.write("\t".join(line) + "\n")
    return trit_file


def export_trit_dir(data, option_name, dir_path, bind_object=None):
    trit_coll = db["sg_feature_file"]
    trit_detail_coll = db["sg_feature_file_detail"]
    trit_dir = os.path.join(dir_path, "trait_dir")
    if not os.path.exists(trit_dir):
        os.mkdir(trit_dir)
    result = trit_coll.find_one({"_id": ObjectId(data)})
    if not result:
        raise Exception("没有在表sg_feature_file中找到id:{}对应的结果".format(data))
    feature_dict = result["feature_dict"]
    specimen_ids = result["specimen_list"]
    task_id = result["task_id"]
    task_result = db["sg_task"].find_one({"task_id": task_id})
    pop_type = task_result["poptype"]
    specimen_list = []
    child_coll = db["sg_specimen_child"]
    results = child_coll.find({"task_id": task_id})
    if not results:
        raise Exception("没有在表sg_specimen_child中找到task_id:{}对应的结果".format(task_id))
    for result in results:
        analysis_name = result["analysis_name"]
        if analysis_name in specimen_ids:
            specimen_list.append(analysis_name)
    specimen_list = list(set(specimen_list))
    trit_dict, trit_list, value_list = {}, [], []
    for trit in feature_dict.keys():
        trit_dict[trit] = []
        trit_list.append(trit)
        value_list.append(feature_dict[trit])
    for sample in specimen_list:
        result = trit_detail_coll.find_one({"feature_file_id": ObjectId(data), "sample": sample})
        if not result:
            raise Exception("没有在表sg_feature_file_detail表中找到样本：{}对应的结果".format(sample))
        sample = result["sample"]
        for t in trit_list:
            if result[t] == "--" or result[t] == "-":
                result[t] = "NaN"
            trit_dict[t].append(str(result[t]))
    for t in trit_list:
        trit_file = os.path.join(trit_dir, feature_dict[t] + ".txt")
        with open(trit_file, "w") as w:
            if pop_type == "F1":
                w.write("ntrt = 1\t\n" + "nind = " + str(len(specimen_list)) + "\t\nmiss = *\t\nsampleID\t" + feature_dict[t] + "\n")
                for i in range(len(specimen_list)):
                    w.write(specimen_list[i] + "\t" + trit_dict[t][i] + "\n")
            else:
                w.write("Genotype," + ",".join(specimen_list) + "\n")
                w.write(feature_dict[t] + "," + ",".join(trit_dict[t]) + "\n")
    return trit_dir


def export_sample_file(data, option_name, dir_path, bind_object=None):
    trit_coll = db["sg_feature_file"]
    trit_detail_coll = db["sg_feature_file_detail"]
    sample_file = os.path.join(dir_path, "sample.txt")
    result = trit_coll.find_one({"_id": ObjectId(data)})
    if not result:
        raise Exception("没有在表sg_feature_file中找到id:{}对应的结果".format(data))
    specimen_ids = result["specimen_list"]
    task_id = result["task_id"]
    specimen_list = []
    # child_coll = db["sg_specimen_child"]
    child_coll = db["sg_specimen"]
    results = child_coll.find({"task_id": task_id})
    if not results:
        raise Exception("没有在表sg_specimen中找到task_id:{}对应的结果".format(task_id))
    for result in results:
        # analysis_name = result["analysis_name"]
        analysis_name = result["initial_name"]
        if analysis_name in specimen_ids:
            specimen_list.append(analysis_name)
    specimen_list = list(set(specimen_list))
    with open(sample_file, "w") as w:
        w.write("\n".join(specimen_list) + "\n")
    return sample_file


def export_nocp_trit_file(data, option_name, dir_path, bind_object=None):
    """
    qingmei
    20180712
    传参形式固定为data, option_name, dir_path, bind_object=None
    controller里传过来sg_feature_file_id为string
    """
    nocp_trit_file = os.path.join(dir_path, "nocp_trit.xls")
    feature_result = db["sg_feature_file"].find_one({"_id": ObjectId(data)})
    feature_samples_result = db['sg_feature_file_detail'].find({"feature_file_id": ObjectId(data)})
    trait_dict = feature_result['feature_dict']   # trait1: trur_name
    lines = {}
    file_header = ['Genotype']
    for record in feature_samples_result:
        if record['sample'] not in file_header:
            file_header.append(record['sample'])
        for trait_name in trait_dict.keys():
            true_trait_name = trait_dict[trait_name]
            if true_trait_name not in lines.keys():
                lines[true_trait_name] = [true_trait_name]  # 每一行的开头为性状名
            trit_value = record[trait_name]
            trit_value = trit_value.strip()
            lines[true_trait_name].append(trit_value)
    # print(lines)
    with open(nocp_trit_file, 'w') as fw:
        fw.write(','.join(file_header) + '\n')
        trit_name = lines.keys()
        trit_name.sort()
        for true_name in trit_name:
            fw.write(','.join(lines[true_name]) + '\n')
    return nocp_trit_file

def export_cp_trit_file(data, option_name, dir_path, bind_object=None):
    """
    qingmei
    20180807
    """
    cp_trit_file = os.path.join(dir_path, "cp_trit.xls")
    feature_result = db["sg_feature_file"].find_one({"_id": ObjectId(data)})
    feature_samples_result = db['sg_feature_file_detail'].find({"feature_file_id": ObjectId(data)})
    trait_dict = feature_result['feature_dict']   # trait1: trur_name
    print(trait_dict)
    lines = {}  # 存储性状名：性状值
    sample_list = []    # 存储样本名
    trit = []   # cp群体上传多个性状值
    for record in feature_samples_result:
        if record['sample'] not in sample_list:
            sample_list.append(record['sample'])
        for trait_name in trait_dict.keys():
            true_trait_name = trait_dict[trait_name]
            if true_trait_name not in trit:
                trit.append(true_trait_name)    # 存储性状值
            if record['sample'] not in lines.keys():
                lines[record['sample']] = [record['sample']]  # 每一行的开头为样本名
            trit_value = record[trait_name]
            trit_value = trit_value.strip()
            lines[record['sample']].append(trit_value)
    with open(cp_trit_file, 'w') as fw:
        num = str(len(trit))
        fw.write("ntrt=" + num + '\n')
        num = str(len(sample_list))
        fw.write("nind=" + num + '\n')
        fw.write("miss=*" + '\n')
        fw.write('sampleID' + '\t' + '\t'.join(trit) + '\n')
        sample_name = lines.keys()
        sample_name.sort()
        for true_name in sample_name:
            fw.write('\t'.join(lines[true_name]) + '\n')
    return cp_trit_file


# data = '5b2c68aaa4e1af2ac8a1db82'   # sg_feature_id
# option_name = ''
# dir_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/"
# export_nocp_trit_file(data, option_name, dir_path)
# data = '5b68ef13c6598d9a688b4568'   # sg_feature_id
# option_name = ''
# dir_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/"
# export_cp_trit_file(data, option_name, dir_path)
# data = "5b2b334fa4e1af24879372a8"
# option_name = None
# dir_path = "/mnt/ilustre/users/sanger-test/test_tmp/hongdong_test"
# export_trit_file(data, option_name, dir_path, bind_object=None)
# export_trit_dir(data="5b722de3c6598d91158b4567", option_name="", dir_path=dir_path, bind_object=None)
