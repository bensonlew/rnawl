# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
import os
import json
from biocluster.config import Config
from bson.objectid import ObjectId

client = Config().get_mongo_client(mtype="meta")
db = client[Config().get_mongo_dbname("meta")]


def export_otu_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组信息(group_detail)获取OTU表
    使用时确保你的workflow的option里group_detail这个字段
    """
    # db = Config().mongo_client[Config().MONGODB]
    table_path = os.path.join(dir_path, "otu_table.xls")
    rep_path = os.path.join(dir_path, "otu_reps.fasta")
    group_path = os.path.join(dir_path, "group_detail.txt")
    bind_obj.logger.debug("正在导出OTU表格文件，路径:%s" % (table_path))
    bind_obj.logger.debug("正在导出OTU fasta文件，路径:%s" % (rep_path))
    my_collection = db['sg_otu_specimen']
    my_results = my_collection.find({"otu_id": ObjectId(data)})
    if not my_results.count():
        raise Exception("otu_id: {}在sg_otu_specimen表中未找到！".format(data))
    samples = list()
    table_dict = {}
    group_detail = bind_obj.sheet.option("group_detail")
    group_method = bind_obj.sheet.option("group_method")
    bind_obj.logger.debug(group_detail)
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            raise Exception(
                "生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        raise Exception(
            "生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    sample_table = db['sg_specimen']
    group_dict = {}
    for k in table_dict:
        group = []
        for sp_id in table_dict[k]:
            sp = sample_table.find_one({"_id": ObjectId(sp_id)})
            if not sp:
                raise Exception(
                    "group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, 'sg_specimen'))
            else:
                samples.append(sp["specimen_name"])
                group.append(sp["specimen_name"])
        group_dict[k] = group
    with open(group_path, "wb") as w:
        w.write(json.dumps(group_dict))
    collection = db['sg_otu_detail']
    with open(table_path, "wb") as f, open(rep_path, "wb") as w:
        header = "OTU ID"
        if group_method == "":
            header += "\t%s" % "\t".join(samples)
        else:
            for g in group_dict:
                header += "\t%s" % g
        f.write(header + "\n")
        for col in collection.find({"otu_id": ObjectId(data)}):
            table_line = "%s\t" % col["otu"]
            if group_method == "":
                for s in samples:
                    table_line += "%s\t" % col[s]
            if group_method == "sum":
                for g in group_dict:
                    otu = 0
                    for s in group_dict[g]:
                        otu += col[s]
                    table_line += "%s\t" % otu
            if group_method == "average":
                for g in group_dict:
                    sample_count = 0
                    otu = 0
                    for s in group_dict[g]:
                        otu += col[s]
                        sample_count += 1
                    table_line += "%s\t" % (float(otu) / sample_count)
            if group_method == "middle":
                for g in group_dict:
                    otu_list = []
                    otu = 0
                    for s in group_dict[g]:
                        otu_list.append(col[s])
                    otu_list = sorted(otu_list)
                    length = len(otu_list)
                    if (length % 2) == 1:
                        otu = otu_list[length / 2]
                    else:
                        otu = float(otu_list[length / 2] +
                                    otu_list[length / 2 - 1]) / 2
                    table_line += "%s\t" % otu
            f.write("%s\n" % table_line)
            line = ">%s\n" % col["otu"]
            line += col["otu_rep"]
            w.write("%s\n" % line)
    paths = ','.join([table_path, rep_path, group_path])
    return paths
