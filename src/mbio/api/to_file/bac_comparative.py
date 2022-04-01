# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import os
import re
import json
from types import StringTypes
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import defaultdict
from mainapp.models.mongo.bac_comparative import BacComparative
from mainapp.controllers.project.bac_comparative_controller import BacComparativeController


def get_mogno():
    client = Config().get_mongo_client(mtype="bac_comparative")
    db = client[Config().get_mongo_dbname("bac_comparative")]
    return db, client

def export_data_from_sample_detail(samples, dir_path, bind_obj=None):
    """
    工作流获取输入文件以及参数
    :param data:
    :return:
    """
    db, client = get_mogno()
    samples = samples.split(",")
    file_path = os.path.join(dir_path, "sample_info.xls")
    task_id = bind_obj['id']
    print(task_id)
    res = db['sample'].find_one({'task_id': task_id})
    print(res)
    main_id = res['_id']
    gff_dir = ''
    if 'gff_status' in res.keys():
        gff_dir = res['gff_path']
    with open(file_path, "w") as f:
        f.write("sample\tspecies\tstrain\tseq_status\tg_location\tgff_file\tgff_path\tseq_file\tseq_path\tstatus\trrna_num\ttrna_num\n")
        for id in samples:
            result = db['sample_detail'].find_one({"sample_id": ObjectId(main_id), "_id": ObjectId(id)})
            print(result)
            if result:
                if 'gff_path' not in result.keys():
                    gff_path = "-"
                else:
                    if re.search(r";",result['gff_file']):
                        gff = []
                        for i in result['gff_file'].split(";"):
                            de = gff_dir + i
                            gff.append(de)
                        gff_path = ";".join(gff)
                    else:
                        bsa = result['gff_file']
                        gff = gff_dir + bsa
                        gff_path = gff
                if 'seq_path' not in result.keys():
                    seq_path = "-"
                else:
                    seq_path = result['seq_path']
                if 'status' not in result.keys():
                    status = "-"
                else:
                    status = result['status']
                f.write("\t".join([result['sample'],result['species'],result['strain'],result['seq_status'],result['g_location'],result['gff_file'],gff_path,result['seq_file'],seq_path,status,result['rrna_no'],result['trna_no']])+"\n")
            else:
                raise Exception("{}在sample_detail中找不到数据！".format(id))
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


def export_group_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detail这个字段
    data 为task_id ; 所有的分组包括all/All/ALL都需要用该函数导出group文件 add by zhujuan 2017.11.09
    """
    db, client = get_mogno()
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    # data = _get_objectid(data)
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
        if data in ["all", "All", "ALL"]:
            f.write("#sample\t" + "##empty_group##" + "\n")
        else:
            f.write("#sample\t" + schema_name + "\n")
    with open(file_path, "ab") as f:
        for k in sorted(table_dict):
            for sp_id in table_dict[k]:
                sample_name = sp_id ## 前端直接传过来是改名后的样本名称
                if not sample_name:
                    raise Exception("group_detail中的样本_id:{}在样本集中未找到".format(sp_id))
                if len(table_dict) == 1 and k in ['all', 'All', 'ALL']:
                    f.write("{}\t{}\n".format(sample_name, sample_name))
                else:
                    f.write("{}\t{}\n".format(sample_name, k))
    client.close()
    return file_path


def export_group_file(table_id, options_name, dir_path, bind_obj=None):
    db, client = get_mogno()
    file_path = os.path.join(dir_path, 'group.txt')
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (options_name, file_path))
    table_id = _get_objectid(table_id)
    group_info = db['specimen_group'].find_one({'_id': table_id})

    category_names = group_info['category_names']
    specimen_names = group_info['specimen_names']

    with open(file_path, 'w') as f:
        f.write("#name\tgroup_name\n")
        for i in range(len(specimen_names)):
            for sp in specimen_names[i]:
                f.write(sp + '\t' + category_names[i] + '\n')
    return file_path
