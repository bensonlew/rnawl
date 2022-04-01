#-*- coding: utf-8 -*-
from web.utils import requeue
from biocluster.config import Config
from mbio.packages.metabolome.relation_config import Config as rl_config
import pandas as pd
from bson.objectid import ObjectId


col_info = {
    # proj_type: {col_type: [main_col, tetail_col, main_id_key, row_name_key]}
    "prok_rna": {
        "exp": ["sg_exp", "sg_exp_detail", "exp_id"],
        "geneset": ["sg_geneset", "sg_geneset_detail", "geneset_id"],
        "kegg_anno": ["sg_annotation_kegg", "sg_annotation_kegg_table", "kegg_id"],
        "gene_diff": ["sg_diff", "sg_diff_detail", "diff_id"],
        "gene_name": ["sg_annotation_query", "sg_annotation_query_detail", "query_id"], # gene_id
        "group_sample": ["sg_specimen_group_compare", "sg_specimen_group"],
        },
    "ref_rna_v2": {
        "exp": ["sg_exp", "sg_exp_detail", "exp_id"],
        "geneset": ["sg_geneset", "sg_geneset_detail", "geneset_id"],
        "kegg_anno": ["sg_annotation_kegg", "sg_annotation_kegg_table", "kegg_id"],
        "gene_diff": ["sg_diff", "sg_diff_detail", "diff_id"],
        "gene_name": ["sg_annotation_query", "sg_annotation_query_detail", "query_id"], # gene_id
        "group_sample": ["sg_specimen_group_compare", "sg_specimen_group"], # gene_id
        },
    "denovo_rna_v2": {
        "exp": ["sg_exp", "sg_exp_detail", "exp_id"],
        "geneset": ["sg_geneset", "sg_geneset_detail", "geneset_id"],
        "kegg_anno": ["sg_annotation_kegg", "sg_annotation_kegg_table", "kegg_id"],
        "gene_diff": ["sg_diff", "sg_diff_detail", "diff_id"],
        "group_sample": ["sg_specimen_group_compare", "sg_specimen_group"],
        },
    "medical_transcriptome": {
        "exp": ["sg_exp", "sg_exp_detail", "exp_id"],
        "geneset": ["sg_geneset", "sg_geneset_detail", "geneset_id"],
        "kegg_anno": ["sg_annotation_kegg", "sg_annotation_kegg_table", "kegg_id"],
        "gene_diff": ["sg_diff", "sg_diff_detail", "diff_id"],
        "gene_name": ["sg_exp", "sg_exp_detail", "exp_id"], # gene_id
        "group_sample": ["sg_specimen_group_compare", "sg_specimen_group"],
        },
    "whole_transcriptome": {
        "exp": ["exp", "exp_detail", "exp_id"],
        "geneset": ["geneset", "geneset_detail", "geneset_id"],
        "kegg_anno": ["annotation_kegg", "annotation_kegg_table", "kegg_id"],
        "gene_diff": ["diff", "diff_detail", "diff_id"],
        "gene_name": ["exp", "exp_detail", "exp_id"],
        "group_sample": ["specimen_group_compare", "specimen_group"],
        },
}


def get_db(proj_type, task_id, db_version=1):
    if db_version:
        cf = Config()
    else:
        cf = rl_config()
    client = cf.get_mongo_client(proj_type, task_id=task_id)
    db = client[cf.get_mongo_dbname(proj_type, task_id=task_id)]
    return db


def dump_trans_exp(detail_col, main_id, main_key):
    ## gene_id 样本，其它信息
    keys = detail_col.find_one({main_key: main_id}).keys()
    if "gene_id" in keys:
        g_id = "gene_id"
    else:
        g_id = "seq_id"
    [keys.remove(k) for k in ["_id", main_key, g_id]]
    keys = [g_id,] + keys
    out_data = []
    for one in detail_col.find({main_key: main_id}).batch_size(1000):
        out_data.append([one[k] if k in one else '-' for k in keys])
    keys[0] = "gene_id"
    df = pd.DataFrame(out_data, columns=keys)
    return df


def dump_trans_geneset(detail_col, main_id, main_key):
    # gene_id regulate(up or down)
    # gene1 up
    # gene2 down
    out_data = []
    for one in detail_col.find({main_key: main_id}).batch_size(1000):
        seq_list = one["seq_list"]
        regulate_list = one["regulate_list"] if "regulate_list" in one else ['-'] * len(seq_list)
        for i in range(len(seq_list)):
            out_data.append([seq_list[i], regulate_list[i]])
    df = pd.DataFrame(out_data, columns=["gene_id", "regulate"])
    return df


def dump_trans_kegg_anno(detail_col, main_id, main_key):
    # gene_id pathways K_lists
    # gene1 map00010;map00011 K01010,K01011
    out_data = []
    filter_d = {main_key: main_id, "seq_type": "all", "anno_type": "G"}
    if not detail_col.find_one(filter_d):
        filter_d = {main_key: main_id, "anno_type": "G"}
    if not detail_col.find_one(filter_d):
        filter_d = {main_key: main_id}
    for one in detail_col.find(filter_d).batch_size(1000):
        out_data.append([one["transcript_id"], one["paths"], one["ko_id"]])
    df = pd.DataFrame(out_data, columns=["gene_id", "pathways", "K_lists"])
    return df


def dump_trans_gene_diff(detail_col, main_id, main_key):
    keys = detail_col.find_one({main_key: main_id}).keys()
    if "gene_id" in keys:
        g_id = "gene_id"
    else:
        g_id = "seq_id"
    [keys.remove(k) for k in ["_id", main_key, g_id]]
    keys = [g_id,] + keys
    out_data = []
    for one in detail_col.find({main_key: main_id}).batch_size(1000):
        out_data.append([one[k] for k in keys])
    keys[0] = "gene_id"
    df = pd.DataFrame(out_data, columns=keys)
    return df

def dump_gene_name(detail_col, main_id, main_key):
    out_data = []
    keys = ["gene_id", "gene_name"]
    for one in detail_col.find({main_key: main_id}).batch_size(1000):
        out_data.append([one[k] if k in one else '' for k in keys])
    df = pd.DataFrame(out_data, columns=keys)
    return df

def dump_group_samples(proj_type, task_id, compare_name, db_version=1):
    '''
    获取差异分析中使用的比较组的样本分组信息
    '''
    db = get_db(proj_type, task_id, db_version=1)
    compare_col, group_col = col_info[proj_type]["group_sample"]
    compare_info = db[compare_col].find_one({"task_id": task_id,
                                             "compare_names": {"$regex": "\"{}\"".format(compare_name)}})
    if not compare_info:
        raise Exception("未找到比较组信息: " + compare_name)
    group_ojb_id = ObjectId(compare_info["specimen_group_id"])

    group_infos = db[group_col].find_one({"$or":[{"main_id": group_ojb_id},
                                                 {"_id": group_ojb_id}]})
    group_names = compare_name.split("|")
    group_samples = {}
    for i in range(len(group_infos["category_names"])):
        if group_infos["category_names"][i] in group_names:
            group_samples[group_infos["category_names"][i]] = group_infos["specimen_names"][i]
    return group_samples

def get_main_id(db, main_col_name, task_id, main_id):
    print("***get_main_id: {} {} {} {}".format(db, main_col_name, task_id, main_id))
    if main_id:
        main_id = ObjectId(main_id)
        main_info = db[main_col_name].find_one({"$or":[{"main_id": main_id}, {"_id": main_id}]})
        print("main_info: {}".format(db[main_col_name]))
        if "batch_main_id" in main_info:
            main_id = main_info["batch_main_id"]
        return main_id, main_info["name"]

    if main_col_name.endswith('exp'):
        filter_d = {"task_id": task_id, "$or":[{"exp_level": 'G'}, {"level": 'G'}]}
    else:
        filter_d = {"task_id": task_id}
    main_info = db[main_col_name].find_one(filter_d)
    if main_info:
        return main_info["_id"], main_info["name"]
    else:
        raise Exception("未找到数据")


dump_funs = {
    "exp": dump_trans_exp,  # 表达数据
    "geneset": dump_trans_geneset,  # 基因集
    "kegg_anno": dump_trans_kegg_anno,  # kegg 注释
    "gene_diff": dump_trans_gene_diff,  # 差异数据表表
    "gene_name": dump_gene_name,  # 基因组名称
    "group_samples": dump_group_samples,  # 分组名称
}


def dump_trans_data(outfile, proj_type, task_id, col_type, main_id=None, compare_name=None, db_version=1):
    if proj_type == "denovo_rna_v2" and col_type == "gene_name":
        return None if outfile else (None, None)
    if col_type == "group_samples":
        dump_group_samples(proj_type, task_id, compare_name)
    db = get_db(proj_type, task_id, db_version)
    main_col_name, detail_col_name, main_id_key = col_info[proj_type][col_type]
    main_id, name = get_main_id(db, main_col_name, task_id, main_id)
    print("***return: {}".format(main_id))
    detail_col = db[detail_col_name]
    df = dump_funs[col_type](detail_col, ObjectId(main_id), main_id_key)
    if outfile:
        df.to_csv(outfile, sep='\t', header=True, index=False)
        return name
    else:
        return df, name
