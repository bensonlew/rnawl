# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
import os
from biocluster.config import Config
from bson.objectid import ObjectId
import types
import json
import re
from types import StringTypes


client = Config().get_mongo_client(mtype="ref_rna")
db = client[Config().get_mongo_dbname("ref_rna")]


def export_express_matrix(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_rna"]
    fpkm_path = os.path.join(dir_path, "%s_fpkm.matrix" % option_name)
    count_path = os.path.join(dir_path, "%s_count.matrix" % option_name)
    bind_obj.logger.debug("正在导出计数矩阵:%s；fpkm矩阵:%s" % (count_path, fpkm_path))
    collection = db['sg_denovo_express_detail']
    my_collection = db['sg_denovo_express']
    results = collection.find({'$and': [{'express_id': ObjectId(data)}, {'type': 'gene'}]})
    my_result = my_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        raise Exception("意外错误，express_id:{}在sg_denovo_express中未找到！".format(ObjectId(data)))
    samples = my_result['specimen']
    with open(fpkm_path, "wb") as f, open(count_path, 'wb') as c:
        head = '\t'.join(samples)
        f.write('\t' + head + '\n')
        c.write('\t' + head + '\n')
        for result in results:
            gene_id = result['gene_id']
            fpkm_write = '{}'.format(gene_id)
            count_write = '{}'.format(gene_id)
            for sam in samples:
                fpkm = sam + '_fpkm'
                count = sam + '_count'
                fpkm_write += '\t{}'.format(result[fpkm])
                count_write += '\t{}'.format(result[count])
            fpkm_write += '\n'
            count_write += '\n'
            f.write(fpkm_write)
            c.write(count_write)
    paths = ','.join([fpkm_path, count_path])
    return paths


def export_control_file(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_rna"]
    file_path = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出计数矩阵:%s" % file_path)
    collection = db['sg_denovo_control']
    result = collection.find_one({'_id': ObjectId(data)})
    if not result:
        raise Exception("意外错误，control_id:{}在sg_denovo_control中未找到！".format(ObjectId(data)))
    group_id = result['group_id']
    if group_id not in ['all', 'All', 'ALL']:
        if isinstance(group_id, types.StringTypes):
            group_id = ObjectId(group_id)
        group_coll = db['sg_denovo_specimen_group']
        g_result = group_coll.find_one({'_id': group_id})
        if not g_result:
            raise Exception("意外错误，control_file的group_id:{}在sg_denovo_specimen_group中未找到！".format(group_id))
    control_detail = result['control_names']
    with open(file_path, 'wb') as w:
        w.write('#control\t{}\n'.format(result['scheme_name']))
        for i in control_detail:
            w.write('{}\t{}\n'.format(i.keys()[0], i.values()[0]))
    return file_path


def export_group_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detal这个字段
    """
    # db = Config().mongo_client[Config().MONGODB + "_rna"]
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    if data in ["all", "All", "ALL"]:
        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
        return file_path
    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')
    group_table = db['sg_denovo_specimen_group']
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        raise Exception("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    group_schema = group_table.find_one({"_id": ObjectId(data)})
    if not group_schema:
        raise Exception("无法根据传入的group_id:{}在sg_denovo_specimen_group表里找到相应的记录".format(data))
    schema_name = re.sub("\s", "_", group_schema["group_name"])  # 将分组方案名的空格替换成下划线
    with open(file_path, "wb") as f:
        f.write("#sample\t" + schema_name + "\n")
    sample_table_name = 'sg_denovo_specimen'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    raise Exception("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen_name"]
                f.write("{}\t{}\n".format(sp_name, k))
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


def export_bam_path(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_rna"]
    my_collection = db['sg_denovo_express']
    my_result = my_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        raise Exception("意外错误，express_id:{}在sg_denovo_express中未找到！".format(ObjectId(data)))
    bam_dir = my_result['bam_path']
    dir_path = bam_dir
    return dir_path


def export_bed_path(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_rna"]
    my_collection = db['sg_denovo_orf']
    my_result = my_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        raise Exception("意外错误，id:{}在sg_denovo_orf中未找到！".format(ObjectId(data)))
    bam_dir = my_result['orf-bed']
    dir_path = bam_dir
    return dir_path


def export_fasta_path(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_rna"]
    my_collection = db['sg_denovo_sequence']
    my_result = my_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        raise Exception("意外错误，id:{}在sg_denovo_sequence中未找到！".format(ObjectId(data)))
    gene_path = my_result['gene_path']
    dir_path = gene_path
    return dir_path


def export_kegg_table(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_rna"]
    kegg_path = os.path.join(dir_path, 'gene_kegg_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg_path))
    par_collection = db["sg_denovo_express_diff"]
    task_id = par_collection.find_one({"_id": ObjectId(data)})["task_id"]
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\n')
        anno_id = db['sg_denovo_annotation'].find_one({'task_id': task_id})['_id']
        results = db['sg_denovo_annotation_kegg_table'].find({'$and': [{'annotation_id': anno_id}, {'type': 'gene'}]})
        if not results:
            raise Exception("生成kegg_table出错：annotation_id:{}在sg_denovo_annotation_kegg_table中未找到！".format(ObjectId(anno_id)))
        for result in results:
            w.write('{}\t{}\t{}\t{}\t{}\n'.format(result['query_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths']))
    return kegg_path


def export_go2level(data, option_name, dir_path, bind_obj=None):
    go2level = os.path.join(dir_path, "go2level.xls")
    bind_obj.logger.debug("正在导出差异基因列表:%s" % go2level)
    par_collection = db["sg_denovo_express_diff"]
    task_id = par_collection.find_one({"_id": ObjectId(data)})["task_id"]
    my_result = db["sg_denovo_annotation"].find_one({"task_id": task_id})
    anno_id = my_result["_id"]
    if not my_result:
        raise Exception("意外错误，annotation_id:{}在sg_denovo_annotation中未找到！".format(anno_id))
    collection = db["sg_denovo_annotation_go_graph"]
    results = collection.find({"$and": [{"annotation_id": anno_id}, {"type": "gene"}, {"level": 2}]})
    if not results:
        raise Exception("生成go2level出错：annotation_id:{}在sg_denovo_annotation_go_graph中未找到！".format(ObjectId(anno_id)))
    with open(go2level, "wb") as w:
        w.write("term_type\tterm\tGO\tnumber\tpercent\tsequence\n")
        for result in results:
            w.write(result['parent_name'] + '\t' + result['go_name'] + '\t' + result['go_id'] + '\t' + str(result['num']) + '\t' + str(result['rate']) + '\t' + result['sequence'] + '\n')
    return go2level

def export_gos_list(data, option_name, dir_path, bind_obj=None):
    gos_list = os.path.join(dir_path, "unigene_GO.list")
    bind_obj.logger.debug("正在导出差异基因go列表:%s" % gos_list)
    par_collection = db["sg_denovo_express_diff"]
    task_id = par_collection.find_one({"_id": ObjectId(data)})["task_id"]
    my_result = db["sg_denovo_annotation"].find_one({"task_id": task_id})
    anno_id = my_result["_id"]
    if not my_result:
        raise Exception("意外错误，annotation_id:{}在sg_denovo_annotation中未找到！".format(anno_id))
    collection = db["sg_denovo_annotation_gos_list"]
    results = collection.find({"annotation_id": ObjectId(anno_id)})
    if not results:
        raise Exception("生成gos_list出错：annotation_id:{}在sg_denovo_annotation_gos_list中未找到！".format(ObjectId(anno_id)))
    with open(gos_list, "wb") as w:
        for result in results:
            gene_id = result["gene_id"]
            go_list = result["go_id"]
            w.write(gene_id + "\t" + go_list + "\n")
    return gos_list

def export_all_gene_list(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_rna"]
    all_list = os.path.join(dir_path, "all_gene.list")
    bind_obj.logger.debug("正在导出所有基因列表:%s" % all_list)
    par_collection = db["sg_denovo_express_diff"]
    express_id = par_collection.find_one({"_id": ObjectId(data)})["express_id"]
    collection = db["sg_denovo_express_detail"]
    my_collection = db["sg_denovo_express"]
    results = collection.find({"$and": [{"express_id": ObjectId(express_id)}, {"type": "gene"}]})
    my_result = my_collection.find_one({"_id": ObjectId(express_id)})
    if not my_result:
        raise Exception("意外错误，expree_id:{}在sg_denovo_express中未找到!".format(ObjectId(express_id)))
    with open(all_list, "wb") as w:
        for result in results:
            gene_id = result["gene_id"]
            w1.write(gene_id + "\n")
        collection1 = db["sg_denovo_express_diff"]
        results1 = collection1.find({"express_id": ObjectId(data)})
        # results1 = collection1.find_one({"express_id": data})
        for result1 in results1:
            express_diff_id = result1["_id"]
            collection2 = db["sg_denovo_express_diff_detail"]
            results2 = collection2.find({"express_diff_id": express_diff_id})
            for result2 in results2:
                gene_id = result2["gene_id"]
                w2.write(gene_id + "\n")
    my_collection1 = db["sg_denovo_annotation"]
    my_result1 = my_collection1.find_one({"task_id": task_id})
    annotation_id = my_result1["_id"]
    if not my_result1:
        raise Exception("意外错误，annotation_id:{}在sg_denovo_annotation中未找到！".format(annotation_id))
    collection3 = db["sg_denovo_annotation_gos_list"]
    results3 = collection3.find({"annotation_id": ObjectId(annotation_id)})
    with open(gos_list, "wb") as w4:
        for result3 in results3:
            gene_id = result3["gene_id"]
            go_list = result3["go_id"]
            w4.write(gene_id + "\t" + go_list + "\n")
    paths = ','.join([all_list, diff_list, gos_list])
    return paths

def export_diff_express(data, option_name, dir_path, bind_obj=None):
    name = bind_obj.sheet.option("name")
    compare_name = bind_obj.sheet.option("compare_name")
    diff_express = os.path.join(dir_path, "%s_vs_%s.diff.exp.xls" % (name, compare_name))
    bind_obj.logger.debug("正在导出差异基因表达量表:%s" % diff_express)
    collection = db["sg_denovo_express_diff_detail"]
    results = collection.find({"$and": [{"express_diff_id": ObjectId(data)}, {"name": name}, {"compare_name": compare_name}]})
    my_collection = db["sg_denovo_express_diff"]
    my_result = my_collection.find_one({"_id": ObjectId(data)})
    if not my_result:
        raise Exception("意外错误，expree_diff_id:{}在sg_denovo_express_diff中未找到!".format(ObjectId(data)))
    with open(diff_express, "wb") as w:
        w.write('gene_id\t{}_count\t{}_count\t{}_fpkm\t{}_fpkm\tlog2fc({}/{})\tpvalue\tfdr\tsignificant\tregulate\n'.format(name, compare_name, name, compare_name, compare_name, name))
        for result in results:
            gene_id = result["gene_id"]
            try:
                name_count = result["{}_count".format(name)]
                compare_count = result["{}_count".format(compare_name)]
                name_fpkm = result["{}_fpkm".format(name)]
                compare_fpkm = result["{}_fpkm".format(compare_name)]
            except:
                name_count = result["{}_mean_count".format(name)]
                compare_count = result["{}_mean_count".format(compare_name)]
                name_fpkm = result["{}_mean_fpkm".format(name)]
                compare_fpkm = result["{}_mean_fpkm".format(compare_name)]
            try:
                significant = result['significant']
                regulate = result['regulate']
                log = result["log2fc({}/{})".format(compare_name, name)]
                pval = result["pvalue"]
                fdr = result["fdr"]
            except:
                significant = result['Significant']
                regulate = result['Regulate']
                log = result["log2FC({}/{})".format(compare_name, name)]
                pval = result["Pvalue"]
                fdr = result["Fdr"]
            w.write(gene_id + '\t' + str(name_count) + '\t' + str(compare_count) + '\t' + str(name_fpkm) + '\t' + str(compare_fpkm) + '\t' + str(log) + '\t' + str(pval) + '\t' + str(fdr) + '\t' + significant + '\t' + regulate + '\n')
    return diff_express
