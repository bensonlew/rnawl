# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# __version__ = 'v1.0'
# __last_modified__ = '20171101'
"""
从数据库中导出一张otu表
"""
from biocluster.config import Config
from bson.objectid import ObjectId
from types import StringTypes
import re

# client = Config().get_mongo_client(mtype="meta")
# db = client[Config().get_mongo_dbname("meta")]


def export_otu_table_by_level(otu_id, level, path):
    client = Config().get_mongo_client(mtype="meta")
    db = client[Config().get_mongo_dbname("meta")]
    collection = db['sg_otu_specimen']
    otu_id = _get_objectid(otu_id)
    results = collection.find({"otu_id": otu_id})
    if not results.count():
        raise Exception("otu_id: {}在sg_otu_specimen表中未找到！".format(otu_id))
    samples = list()
    for result in results:
        if "specimen_id" not in result:
            raise Exception("otu_id:{}错误，请使用新导入的OTU表的id".format(otu_id))
        sp_id = result['specimen_id']
        my_collection = db['sg_specimen']
        my_result = my_collection.find_one({"_id": sp_id})
        if not my_result:
            raise Exception("样本id:{}在sg_specimen表里未找到".format(sp_id))
        samples.append(my_result["specimen_name"])
    level = int(level)
    collection = db['sg_otu_detail']
    name_dic = dict()
    results = collection.find({"otu_id": otu_id})
    if not results.count():
        raise Exception("otu_id: {}在sg_otu_detail表中未找到！".format(otu_id))
    for col in results:
        tmp = level + 1
        new_classify_name = _create_classify_name(col, tmp)
        if new_classify_name not in name_dic:
            name_dic[new_classify_name] = dict()
            for sp in samples:
                name_dic[new_classify_name][sp] = int(col[sp])
        else:
            for sp in samples:
                name_dic[new_classify_name][sp] += int(col[sp])
    with open(path, "wb") as f:
        f.write("OTU ID\t%s\n" % "\t".join(samples))
        for k in name_dic.iterkeys():
            line = k
            for s in samples:
                line += "\t" + str(name_dic[k][s])
            line += "\n"
            f.write(line)


def _get_objectid(otu_id):
    if not isinstance(otu_id, ObjectId):
        if not isinstance(otu_id, StringTypes):
            raise Exception("{}不为ObjectId类型或者其对应的字符串".format(otu_id))
        else:
            try:
                otu_id = ObjectId(otu_id)
            except:
                raise Exception("{}不为ObjectId类型或者其对应的字符串".format(otu_id))
    return otu_id


def _create_classify_name(col, tmp):
    """
    在数据库读取OTU表的分类信息，形成OTU表的第一列
    """
    LEVEL = {
        1: "d__", 2: "k__", 3: "p__", 4: "c__", 5: "o__",
        6: "f__", 7: "g__", 8: "s__", 9: "otu"
    }
    my_list = list()
    last_classify = ""
    for i in range(1, 10):
        if LEVEL[i] in col:
            if re.search("uncultured$", col[LEVEL[i]]) or re.search("Incertae_Sedis$", col[LEVEL[i]]) or re.search(
                    "norank$", col[LEVEL[i]]):
                if i == 1:
                    pass
                    #raise Exception("在域水平上的分类为uncultured或Incertae_Sedis或是norank")
                else:
                    col[LEVEL[i]] = col[LEVEL[i]] + "_" + col[LEVEL[i - 1]]
    for i in range(1, tmp):
        if LEVEL[i] not in col:
            if last_classify == "":
                last_classify = col[LEVEL[i - 1]]
            my_str = LEVEL[i] + "Unclasified_" + last_classify
        else:
            if not col[LEVEL[i]]:
                if LEVEL[i] == "otu":
                    my_str = "otu__" + "Unclasified_" + last_classify
                else:
                    my_str = LEVEL[i] + "Unclasified_" + last_classify
            else:
                my_str = col[LEVEL[i]]
        my_list.append(my_str)
    new_classify_name = "; ".join(my_list)
    return new_classify_name
