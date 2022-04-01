# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last_modified by guhaidong 20171116
import os
from biocluster.config import Config
from bson.objectid import ObjectId
import json

client = Config().get_mongo_client(mtype="meta")
db = client[Config().get_mongo_dbname("meta")]


def export_est_table(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB]
    est_path = os.path.join(dir_path, "%s_input.estimators.xls" % option_name)
    # file_path = os.path.join(dir_path, "%s_input.est_for_t.xls" % option_name)
    # cmd_path = os.path.join(dir_path, "cmd.r")
    bind_obj.logger.debug("正在导出参数%s的多样性指数表格为文件，路径:%s" % (option_name, est_path))
    collection = db['sg_alpha_diversity_detail']
    est_collection = db['sg_alpha_diversity']
    task_id = "_".join(bind_obj.sheet.id.split('_')[0:2])  # add task_id by guhaidong 20171116
    result = est_collection.find_one({"_id": ObjectId(data), "task_id": task_id})  # add task_id by guhaidong 20171116
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
    # print write_data
    with open(est_path, "wb") as f:
        for col in details:
            f.write("\t{}".format(col["specimen_name"]))
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
    # print("hhhhhhhhhhhh")
    return est_path
