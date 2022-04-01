# -*- coding: utf-8 -*-
from biocluster.config import Config
import pandas as pd
from bson.objectid import ObjectId
import json

### tool 中 import 使用
### 通过主表id从mongo导出 **代谢** 的相关数据


col_info = {
    "anno_hmdb": ["anno_hmdb_level", "hmdb_id"],
    "anno_keggc": ["anno_keggc_level", "kegg_id"],
    "anno_keggp": ["anno_keggp_level", "kegg_id"],
    "anno_overview": ["anno_overview_detail", "overview_id"],
    "exp": ["exp_mix,exp_pos,exp_neg", "exp_id"],
    "exp_diff": ["exp_diff_detail", "diff_id"]
}


def get_db(ref=False):
    cf = Config()
    client = cf.get_mongo_client("metabolome", ref=ref)
    db = client[cf.get_mongo_dbname("metabolome", ref=ref)]
    return db


def dump_mongo_data(outfile, obj_id, col_name, detail_col=None, main_key=None, index_key=None):
    db = get_db()
    if not detail_col:
        detail_col, main_k = col_info[col_name]
    if not main_key:
        main_key = main_k
    df = get_mongo_data(db[detail_col], obj_id, main_key, index_key)
    if outfile:
        df.to_csv(outfile, sep='\t', header=True, index=False)
    else:
        return df


def get_mongo_data(col, obj_id, main_key, index_key=None):
    keys = col.find_one({main_key: ObjectId(obj_id)}).keys()
    keys.remove("_id")
    keys.remove(main_key)
    if index_key:
        keys.remove(index_key)
        keys = [index_key] + keys
    out_data = []
    for one in col.find({main_key: ObjectId(obj_id)}).batch_size(1000):
        out_data.append([one[k] if k in one else '' for k in keys])
    df = pd.DataFrame(out_data, columns=keys)
    return df

def expdiff_group_samples(diff_id):
    db = get_db()
    diff_info = db["exp_diff"].find_one({"_id": ObjectId(diff_id)})

    params = json.loads(diff_info["params"])
    return params["group_detail"]

def get_samples_map(task_id, relation_task_id):
    db = get_db()
    relation_info = db["sg_relation_analysis"].find_one({"task_id": task_id,
                                                         "relate_task_id": relation_task_id})
    if not relation_info:
        print("没找到样本对应关系db: {}, metab: {}, trans: {}".format(db, task_id, relation_task_id))
        raise Exception("没找到样本对应关系 metab: {}, trans: {}".format(task_id, relation_task_id))

    return dict(zip(relation_info["relate_sp_name"], relation_info["sp_name"]))