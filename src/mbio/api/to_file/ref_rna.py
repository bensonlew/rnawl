# -*- coding: utf-8 -*-
# __author__ = 'sanger'
from __future__ import division
import os
from biocluster.config import Config
from bson.objectid import ObjectId
import types
import json
import re
from types import StringTypes
import gridfs
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer

# db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]


def export_gene_list(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    gene_list_path = os.path.join(dir_path, "%s_gene.list" % option_name)
    bind_obj.logger.debug("正在导出基因集")
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(data)))
    results = collection.find_one({"geneset_id": ObjectId(data)})
    with open(gene_list_path, "wb") as f:
        gene_list = results["gene_list"]
        for gene_id in gene_list:
            f.write(gene_id + "\n")
    return gene_list_path


def export_blast_table(data, option_name, dir_path, bind_obj=None):
    """
    获取blast结果
    """
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    nr_table_path = os.path.join(dir_path, "nr_{}.xls".format(option_name))
    gene_nr_table_path = os.path.join(dir_path, "gene_nr_{}.xls".format(option_name))
    sw_table_path = os.path.join(dir_path, "swissprot_{}.xls".format(option_name))
    gene_sw_table_path = os.path.join(dir_path, "gene_swissprot_{}.xls".format(option_name))
    stat_collection = db["sg_annotation_stat"]
    stat_result = stat_collection.find({"_id": ObjectId(data)})
    if not stat_result.count():
        bind_obj.set_error("stat_id:{}在sg_annotation_stat表中未找到".format(data))
    for result in stat_result:
        task_id = result["task_id"]
    blast_collection = db["sg_annotation_blast"]
    blast_result = blast_collection.find({"task_id": task_id, "stat_id": ObjectId(data)})
    if not blast_result.count():
        bind_obj.set_error("stat_id:{}在sg_annotation_blast表中未找到".format(data))
    for result in blast_result:
        blast_id = result["_id"]
    collection = db["sg_annotation_blast_detail"]
    results = collection.find({"blast_id": blast_id})
    with open(nr_table_path, "wb") as w1, open(gene_nr_table_path, "wb") as w2, open(sw_table_path, "wb") as w3, open(gene_sw_table_path, "wb") as w4:
        header = "Score\tE-Value\tHSP-Len\tIdentity-%\tSimilarity-%\tQuery-Name\tQ-Len\tQ-Begin\t"
        header += "Q-End\tQ-Frame\tHit-Name\tHit-Len\tHsp-Begin\tHsp-End\tHsp-Frame\tHit-Description\n"
        w1.write(header)
        w2.write(header)
        w3.write(header)
        w4.write(header)
        for result in results:
            db = result["database"]
            anno_type = result["anno_type"]
            seq_type = result["seq_type"]
            score = result["score"]
            evalue = result["e_value"]
            hsp_len = result["hsp_len"]
            identity = result["identity_rate"]
            similarity = result["similarity_rate"]
            query_id = result["query_id"]
            hit_name = result["hit_name"]
            description = result["description"]
            q_len = result["q_len"]
            q_begin = result["q_begin"]
            q_end = result["q_end"]
            q_frame = result["q_frame"]
            hit_len = result["hit_len"]
            hsp_begin = result["hsp_begin"]
            hsp_end = result["hsp_end"]
            hsp_frame = result["hsp_frame"]
            line = str(score) + "\t" + str(evalue) + "\t" + str(hsp_len) + "\t" + str(identity) + "\t"
            line += str(similarity) + "\t" + query_id + "\t" + str(q_len) + "\t" + q_begin + "\t" + q_end + "\t" + q_frame + "\t"
            line += hit_name + "\t" + str(hit_len) + "\t" + hsp_begin + "\t" + hsp_end + "\t" + hsp_frame + "\t" + description + "\n"
            if seq_type == "new":
                if db == "nr":
                    if anno_type == "transcript":
                        w1.write(line)
                    if anno_type == "gene":
                        w2.write(line)
                if db == "swissprot":
                    if anno_type == "transcript":
                        w3.write(line)
                    if anno_type == "gene":
                        w4.write(line)
    paths = ",".join([nr_table_path, gene_nr_table_path, sw_table_path, gene_sw_table_path])
    return paths


def export_go_list(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    go_list_path = os.path.join(dir_path, "GO.list")
    bind_obj.logger.debug("正在导出%sgo列表:%s" % (option_name, go_list_path))
    geneset_collection = db["sg_geneset"]
    task_id = geneset_collection.find_one({"_id": ObjectId(data)})["task_id"]
    my_result = db["sg_annotation_go"].find_one({"task_id": task_id})
    geneset_type= bind_obj.sheet.option('geneset_type')
    go_id = my_result["_id"]
    if not my_result:
        bind_obj.set_error("意外错误，annotation_go_id:{}在sg_annotation_go中未找到！".format(go_id))
    collection = db["sg_annotation_go_list"]
    results = collection.find({"go_id": ObjectId(go_id)})
    one_record = collection.find_one({"go_id": ObjectId(go_id)})
    if not one_record:
        bind_obj.set_error("生成gos_list出错：annotation_id:{}在sg_annotation_gos_list中未找到！".format(ObjectId(go_id)))
    with open(go_list_path, "wb") as w:
        for result in results:
            gene_id = result["gene_id"]
            go_list = result["gos_list"]
            go_anno_type = result["anno_type"]
            if go_anno_type == geneset_type:
                w.write(gene_id + "\t" + go_list + "\n")
    return go_list_path


def export_kegg_table(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    kegg_path = os.path.join(dir_path, 'gene_kegg_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg_path))
    geneset_collection = db["sg_geneset"]
    bind_obj.logger.debug(data)
    geneset_result = geneset_collection.find_one({"_id": ObjectId(data)})
    task_id = geneset_result["task_id"]
    bind_obj.logger.debug("ttttttt")
    bind_obj.logger.debug(task_id)
    geneset_type = geneset_result["type"]
    # my_result = db["sg_annotation_kegg"].find_one({"task_id": task_id, "seq_type": "new"})
    my_result = db["sg_annotation_kegg"].find({"task_id": task_id})
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\n')
        for main_table in my_result:
            kegg_id = main_table["_id"]
            bind_obj.logger.debug(kegg_id)
            if not my_result:
                bind_obj.set_error("意外错误，annotation_kegg_id:{}在sg_annotation_kegg中未找到！".format(kegg_id))
            # with open(kegg_path, 'wb') as w:
            #     w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\n')
            results = db['sg_annotation_kegg_table'].find({'kegg_id': kegg_id, 'anno_type': geneset_type})
            one_record = db['sg_annotation_kegg_table'].find_one({'kegg_id': kegg_id, 'anno_type': geneset_type})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:{}在sg_annotation_kegg_table中未找到！".format(ObjectId(kegg_id)))
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['ko_id'] + result['query_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'

                w.write('{}\t{}\t{}\t{}\t{}\n'.format(result['query_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths']))
    return kegg_path

def export_kegg_pdf(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    fs = gridfs.GridFS(db)
    annotation_collection = db["sg_annotation_kegg"]
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    main_id = annotation_collection.find_one({"task_id":task_id})["_id"]
    kegg_level_collection = db["sg_annotation_kegg_level"]
    results = kegg_level_collection.find({"kegg_id":main_id, "seq_type":"all", "anno_type":anno_type})
    anno_path = dir_path + "/png"
    if not os.path.exists(anno_path):
        os.mkdir(anno_path)
    for result in results:
        graph_id = result["graph_png_id"]
        pathway_id = result["pathway_id"]
        bind_obj.logger.info("正在导出{}的png文件".format(pathway_id))
        with open(anno_path + "/" + pathway_id + ".png", "w") as fw:
            content = fs.get(graph_id).read()
            fw.write(content)
    return anno_path

def export_all_list(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    all_list = os.path.join(dir_path, "all_gene.list")
    bind_obj.logger.debug("正在导出所有背景基因{}".format(all_list))
    # collection = conn['sg_geneset_detail']
    collection = db['sg_express_detail']
    exp_collection = db['sg_express']
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    task_id = my_result["task_id"]
    geneset_type = my_result["type"]
    bind_obj.logger.debug(task_id)
    # exp_result = exp_collection.find_one({'task_id': task_id, "genes": True, "trans": True})
    exp_result = exp_collection.find_one({'task_id': task_id, "name": {"$regex":"fpkm"}, "trans": True})
    # 获取到rsem的表达量主表
    # my_result = main_collection.find_one({'_id': ObjectId("591aefeba4e1af3ec14249c8")})
    # print exp_result["_id"]
    #############################
    if not exp_result:
        bind_obj.set_error("意外错误，task_id:{}的背景基因在sg_geneset中未找到！".format(data))
    exp_id = exp_result["_id"]
    results = collection.find({"express_id": ObjectId(exp_id), "type": geneset_type, "sample_group": "sample"})
    with open(all_list, "wb") as f:
        for result in results:
            gene_id = result['seq_id']
            f.write(gene_id + "\n")
    return all_list


def export_diff_express(data, option_name, dir_path, bind_obj=None):  # 需要修改
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    name = bind_obj.sheet.option("name")
    compare_name = bind_obj.sheet.option("compare_name")
    diff_express = os.path.join(dir_path, "%s_vs_%s.diff.exp.xls" % (name, compare_name))
    bind_obj.logger.debug("正在导出差异基因表达量表:%s" % diff_express)
    collection = db["sg_denovo_express_diff_detail"]
    results = collection.find({"$and": [{"express_diff_id": ObjectId(data)}, {"name": name}, {"compare_name": compare_name}]})
    my_collection = db["sg_denovo_express_diff"]
    my_result = my_collection.find_one({"_id": ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，expree_diff_id:{}在sg_denovo_express_diff中未找到!".format(ObjectId(data)))
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
            w.write(gene_id + '\t' + str(name_count) + '\t' + str(compare_count) + '\t' + str(name_fpkm) + '\t' + str(compare_fpkm) + '\t' + str(log) + '\t' + str(pval) + '\t' + str(fdr) +
'\t' + significant + '\t' + regulate + '\n')
    return diff_express


def export_cog_class(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    cog_path = os.path.join(dir_path, 'cog_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(cog_path))
    genesets, table_title, task_id, geneset_type = get_geneset_detail(data)
    cog_collection = db["sg_annotation_cog"]
    cog_detail_collection = db["sg_annotation_cog_detail"]
    cog_id = cog_collection.find_one({"task_id": task_id})["_id"]
    # cog_results = cog_detail_collection.find({'cog_id': cog_id})
    cog_results = cog_detail_collection.find({'cog_id': cog_id, 'seq_type': 'all', 'anno_type': geneset_type})
    new_table_title = []
    for tt in table_title:
        # new_tt = [tt + "_COG", tt + "_NOG", tt + "_KOG", tt + "_COG_list", tt + "_NOG_list", tt + "_KOG_list"]
        new_tt = [tt + "_COG", tt + "_NOG", tt + "_COG_list", tt + "_NOG_list"]
        new_table_title = new_table_title + new_tt
    bind_obj.logger.debug(table_title)
    with open(cog_path, "wb") as w:
        w.write("Type\tFunctional Categoris\t" + "\t".join(new_table_title) + "\n")
        for cr in cog_results:
            kog_list = set([])
            nog_list = set(cr["nog_list"].split(";") if cr["nog_list"] else [])
            cog_list = set(cr["cog_list"].split(";") if cr["cog_list"] else [])
            write_line = {}
            for gt in genesets:
                # kog_count = list(kog_list & genesets[gt][1])
                nog_count = list(nog_list & genesets[gt][1])
                cog_count = list(cog_list & genesets[gt][1])
                # if not len(kog_count) + len(nog_count) + len(cog_count) == 0:
                #     write_line[gt] = [str(len(cog_count)), str(len(nog_count)), str(len(kog_count)), ";".join(cog_count), ";".join(nog_count), ";".join(kog_count)]
                if not len(nog_count) + len(cog_count) == 0:
                    write_line[gt] = [str(len(cog_count)), str(len(nog_count)), ";".join(cog_count), ";".join(nog_count)]
            if len(write_line) > 0:
                w.write("{}\t{}\t".format(cr["type"], cr["function_categories"]))
                for tt in table_title:
                    w.write("\t".join(write_line[tt]) + "\t") if tt in write_line else w.write("0\t0\tnone\tnone\t")
                w.write("\n")
    return cog_path


def get_geneset_detail(data):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    geneset_collection = db["sg_geneset"]
    genesets = {}
    names = []
    print data.split(",")
    task_id = ""
    geneset_type = "gene"
    for geneset_id in data.split(","):
        geneset_result = geneset_collection.find_one({"_id": ObjectId(geneset_id)})
        if not geneset_result:
            bind_obj.set_error("意外错误:未找到基因集_id为{}的基因集信息".format(geneset_id))
        task_id = geneset_result["task_id"]
        geneset_type = geneset_result["type"]
        geneset_name = geneset_result["name"]
        names.append(geneset_name)
        genesets[geneset_name] = [geneset_type]
        collection = db['sg_geneset_detail']
        results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
        geneset_names = set(results["gene_list"])
        genesets[geneset_name].append(geneset_names)
    return genesets, names, task_id, geneset_type


def export_go_class(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    go_path = os.path.join(dir_path, 'go_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(go_path))
    genesets, names, task_id, seq_type = get_geneset_detail(data)
    go_collection = db["sg_annotation_go"]
    # go_level_collection = db["sg_annotation_go_level"]
    go_level_collection = db["sg_annotation_go_detail"]
    go_id = go_collection.find_one({"task_id": task_id})["_id"]
    # go_results = go_level_collection.find({'go_id': go_id, "level": 2, "seq_type": "all", "anno_type": seq_type})
    one_record = go_level_collection.find_one({'go_id': go_id, "level": 2, "seq_type": "all", "anno_type": seq_type})
    if not one_record:
        bind_obj.set_error("意外错误:未找到go_id为{}的基因集信息".format(go_id))
    new_table_title = []
    for gt in genesets:
        new_table_title.append(gt + " num")
        new_table_title.append(gt + " percent")
        new_table_title.append(gt + " list")
    bind_obj.logger.debug(new_table_title)
    with open(go_path, "wb") as w:
        w.write("Term type\tTerm\tGO\t" + "\t".join(new_table_title) + "\n")
        term_list = ["biological_process", "cellular_component", "molecular_function"]
        for item in term_list:
            go_results = go_level_collection.find({'go_id': go_id, "level": 2, "seq_type": "all", "anno_type": seq_type})
            for gr in go_results:
                if gr["goterm"] == item:
                    seq_list = set(gr["seq_list"].split(";"))
                    write_line = {}
                    for gt in genesets:
                        total_gene_num = len(genesets[gt][1])
                        go_count = list(seq_list & genesets[gt][1])
                        if not len(go_count) == 0:
                            write_line[gt] = str(len(go_count)) + "\t" + str(len(go_count)/total_gene_num) + "(" + str(len(go_count)) + "/" + str(total_gene_num) + ")" + "\t" + ";".join(go_count)
                    if len(write_line):
                        w.write("{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"]))
                        for tt in genesets:
                            w.write(write_line[tt] + "\t") if tt in write_line else w.write("0\t0\tnone\t")
                        w.write("\n")
    return go_path


def export_gene_list_ppi(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    gene_list_path = os.path.join(dir_path, "%s_list.txt" % option_name)
    bind_obj.logger.debug("正在导出基因集")
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(data)))
    results = collection.find_one({"geneset_id": ObjectId(data)})["gene_list"]
    with open(gene_list_path, "wb") as f:
        f.write("gene_id" + "\n")
        for result in results:
            f.write(result + "\n")
    bind_obj.logger.debug("基因集导出成功！")
    return gene_list_path


# ############表达量部分
####################################################表达量部分
def export_express_matrix_level(data,option_name,dir_path,bind_obj=None):
    """
    type对应的是gene/transcript字段，workflow里确保有这个字段
    express_level对应的是fpkm/tpm字段，workflow里确保有这个字段
    """
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    fpkm_path = os.path.join(dir_path, "%s_fpkm.matrix" % option_name)
    count_path = os.path.join(dir_path, "%s_count.matrix" % option_name)
    bind_obj.logger.debug("正在导出计数矩阵:%s；fpkm矩阵:%s" % (count_path, fpkm_path))
    collection = db['sg_express_detail']
    my_collection = db['sg_express']
    type = bind_obj.sheet.option("type")
    bind_obj.logger.debug(type)
    level = bind_obj.sheet.option("express_level")

    group_detail = bind_obj.sheet.option('group_detail')
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            bind_obj.set_error("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        bind_obj.set_error("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))

    sample_table_name = 'sg_specimen'
    sample_table = db[sample_table_name]
    samples = []
    for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    bind_obj.set_error("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen_name"]
                    samples.append(sp_name)

    #sample_group = bind_obj.sheet.option("sample_group")
    results = collection.find({'$and': [{'express_id': ObjectId(data)}, {'type': '{}'.format(type)},{"sample_group":"sample"},{"value_type":level}]})
    count_results = collection.find({'$and': [{'express_id': ObjectId(data)}, {'type': '{}'.format(type)},{"sample_group":"sample"},{"value_type":"count"}]})
    my_result = my_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，express_id:{}在sg_express中未找到！".format(ObjectId(data)))
    # samples = my_result['specimen']
    def write_file(path, collcetion_results):
        with open(path, "wb") as f:
            head = '\t'.join(samples)
            f.write('\t' + head + '\n')
            for result in collcetion_results:
                #bind_obj.logger.debug(result)
                gene_id = result['seq_id']
                fpkm_write = '{}'.format(gene_id)
                for sam in samples:
                    fpkm = sam
                    try:
                        fpkm_write += '\t{}'.format(result[fpkm])
                        #count_write += '\t{}'.format(result[count])
                    except Exception:
                        print fpkm_write
                        print sam
                        print result
                        #bind_obj.set_error("{}错误".format(result[fpkm]))
                fpkm_write += '\n'
                #count_write += '\n'
                f.write(fpkm_write)
                #c.write(count_write)
    write_file(fpkm_path, results)
    write_file(count_path, count_results)
    paths = ','.join([fpkm_path, count_path])
    return paths

def export_group_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detal这个字段
    """
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    if data in ["all", "All", "ALL"]:
        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
        return file_path
    data = _get_objectid(data, bind_obj=bind_obj)
    group_detail = bind_obj.sheet.option('group_detail')  #另传字段
    group_table = db['sg_specimen_group']
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            bind_obj.set_error("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    if not isinstance(table_dict, dict):
        bind_obj.set_error("生成group表失败，传入的{}不是一个字典或者是字典对应的字符串".format(option_name))
    group_schema = group_table.find_one({"_id": ObjectId(data)})
    print group_schema
    print data
    if not group_schema:
        bind_obj.set_error("无法根据传入的group_id:{}在sg_specimen_group_compare表里找到相应的记录".format(data))
    schema_name = re.sub("\s", "_", group_schema["group_name"])  # 将分组方案名的空格替换成下划线
    with open(file_path, "wb") as f:
        f.write("#sample\t" + "group" + "\n")
    sample_table_name = 'sg_specimen'
    sample_table = db[sample_table_name]
    ## modified by shicaiping at 20180904, purpose to get a sorted group name
    string_k = list()
    for k in table_dict:
        if isinstance(k, unicode):
            string_k.append(k.encode('unicode-escape').decode('string_escape'))
        else:
            string_k.append(k)
    with open(file_path, "ab") as f:
        for k in sorted(string_k, key=str.lower):
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"_id": ObjectId(sp_id)})
                if not sp:
                    bind_obj.set_error("group_detal中的样本_id:{}在样本表{}中未找到".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen_name"]
                f.write("{}\t{}\n".format(sp_name, k))
    return file_path

def _get_objectid(data, bind_obj=None):
    if not isinstance(data, ObjectId):
        if not isinstance(data, StringTypes):
            bind_obj.set_error("{}不为ObjectId类型或者其对应的字符串".format(data))
        else:
            try:
                data = ObjectId(data)
            except:
                bind_obj.set_error("{}不为ObjectId类型或者其对应的字符串".format(data))
    return data

def export_control_file(data, option_name, dir_path, bind_obj=None):  #此函数待定 不一定对
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    file_path = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出计数矩阵:%s" % file_path)
    collection = db['sg_specimen_group_compare']
    result = collection.find_one({'_id': ObjectId(data)})
    if not result:
        bind_obj.set_error("意外错误，control_id:{}在sg_specimen_group_compare中未找到！".format(ObjectId(data)))
    group_id = result['specimen_group_id']
    if group_id not in ['all', 'All', 'ALL']:
        """检查group_id的信息"""
        if isinstance(group_id, types.StringTypes):
            group_id = ObjectId(group_id)
        group_coll = db['sg_specimen_group']
        g_result = group_coll.find_one({'_id': group_id})
        if not g_result:
            bind_obj.set_error("意外错误，control_file的group_id:{}在sg_specimen_group中未找到！".format(group_id))
    control_detail = json.loads(result['compare_names'])
    with open(file_path, 'wb') as w:
        w.write('#control\t{}\n'.format('group'))
        for i in control_detail:    #此处需要修改, 可能会有错误
            # w.write('{}\t{}\n'.format(i.keys()[0], i.values()[0]))
            print i
            control_other = i.split("|")
            w.write('{}\t{}\n'.format(control_other[0], control_other[1]))
    return file_path

def _get_gene_id(geneset,geneset_detail,_id, bind_obj=None):
    try:
        results = geneset_detail.find_one({"geneset_id":ObjectId(_id)})
        seq_id = results['gene_list']
    except Exception:
        bind_obj.set_error("{}在sg_geneset_detail表中没有找到!")
    try:
        my_result = geneset.find_one({"_id":ObjectId(_id)})
        _name = my_result['name']
    except Exception:
        bind_obj.set_error("{}在sg_geneset表中没有找到!")
    return seq_id, _name

def export_geneset_venn_level(data, option_name, dir_path, bind_obj=None):
    """
    level对应的是gene/transcript字段，workflow里确保有这个字段
    """
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    geneset_venn = os.path.join(dir_path,"%s_geneset_venn" %(option_name))
    bind_obj.logger.debug("正在导出计数矩阵:%s" %(geneset_venn))
    collection = db["sg_geneset_detail"]
    my_collection = db["sg_geneset"]
    level = bind_obj.sheet.option("type")
    geneset_table = open(geneset_venn,'w+')
    if re.search(',',data):
        new_geneset_id = data.split(",")
    else:
        new_geneset_id = data
    for ll in new_geneset_id:
        seq,_name = _get_gene_id(geneset=my_collection,geneset_detail=collection,_id = ll)
        _seq = ",".join(seq)
        geneset_table.write(_name+"\t"+_seq+"\n")
    geneset_table.close()
    return geneset_venn

def export_class_code(data,option_name,dir_path,bind_obj=None): #输出class_code信息
    """
    type: 对应的是gene 或transcript
    """
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    class_code = os.path.join(dir_path, "%s_class_code" % option_name)
    bind_obj.logger.debug("正在导出class_code信息:%s" %(class_code))
    type= bind_obj.sheet.option('type')
    class_code_type =  bind_obj.sheet.option("class_code_type")
    class_code_detail = db['sg_express_class_code_detail']
    class_code_col = db['sg_express_class_code']
    class_code_main_info = class_code_col.find_one({"_id":ObjectId(data)})
    task_id = class_code_main_info["task_id"]
    sg_task = db["sg_task"]
    sg_task_info = sg_task.find_one({"task_id":task_id})
    if sg_task_info["is_demo"] == 2:
        demo_id = sg_task_info["demo_id"]  # sanger_21455
        class_code_demo_info = class_code_col.find_one({"task_id":demo_id})
        bind_obj.logger.info(class_code_demo_info)
        data = str(class_code_demo_info["_id"])
    bind_obj.logger.info(data)
    class_code_info = class_code_detail.find({"class_code_id":ObjectId(data),"type":class_code_type})
    with open(class_code,'w+') as f:
        header = ['seq_id','gene_name',"class_code"]
        f.write("\t".join(header)+"\n")
        for d in class_code_info:
            if type == 'gene':
                _write = d['assembly_gene_id']+"\t"+d['gene_name']+"\t" + d["class_code"] + "\n"
            if type == 'transcript':
                _write = d['assembly_trans_id']+"\t"+d['gene_name']+"\t" + d["class_code"] + "\t" + d["assembly_gene_id"] + "\n"
            f.write(_write)
    return class_code

def export_add_info(data,option_name,dir_path,bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出add_info信息")
    col = db["sg_annotation_kegg"]
    result = col.find_one({"task_id":task_id})
    insert_id = result["_id"]
    print insert_id
    print anno_type
    col = db["sg_annotation_kegg_level"]
    results = col.find({"kegg_id":insert_id, "seq_type":"all", "anno_type":anno_type})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info

if __name__ == "__main__":
    data = "5909a269a4e1af11112543e2"
    option_name = "class_code"
    dir_path = "/mnt/ilustre/users/sanger-dev/workspace/20170505/DiffExpress_tsg_1000_4773_2935"
    export_class_code(data,option_name,dir_path)
    print 'end!'

def export_geneset_cluster_level(data,option_name,dir_path,bind_obj=None):  #这个函数待定 并且导表函数也待定
    """
    此函数暂时没有用到
    log对应的是2/10字段，workflow里确保有这个字段
    type对应的是fpkm/tpm字段，workflow里确保有这个字段
    data 是两个id，由逗号连接，第一个id是geneset_id 第二个是express_id
    """
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    fpkm_path = os.path.join(dir_path, "%s_fpkm.matrix" % option_name)
    bind_obj.logger.debug("正在导出表达量矩阵矩阵:%s" %(fpkm_path))
    log = bind_obj.sheet.option("log")
    type = bind_obj.sheet.option('type')
    if not re.search(',',data):
        bind_obj.set_error("{}必须是两个ObjectId对象,并由逗号连接".format(data))
    geneset_id = data.split(",")[0]
    express_id = data.split(",")[1]
    geneset_collection = db['sg_geneset']
    geneset_detail_collection = db['sg_geneset_detail']
    express_collection = db['sg_express']
    express_detail_collection = db['sg_express_detail']
    seq,_name = _get_gene_id(geneset_collection,geneset_detail_collection,geneset_id)
    express_data = express_collection.find_one({"_id":ObjectId(express_id)})
    samples = express_data["specimen"]
    with open(fpkm_path,"w+") as f:
        head = '\t'.join(samples)
        f.write('\t' + head + '\n')
        for seq_id in seq:
            out = express_detail_collection.find_one({'$and':[{"express_id":ObjectId(express_id),"seq_id":seq_id}]})
            print out
            fpkm_write = '{}'.format(seq_id)
            for sam in samples:
                if log ==2:
                    fpkm = sam + '_log2_{}'.format(type)
                elif log == 10:
                    fpkm = sam + '_log10_{}'.format(type)
                else:
                    fpkm = sam + '_{}'.format(type)
                print fpkm
                print out[fpkm]
                fpkm_write += '\t{}'.format(out[fpkm])
            fpkm_write += '\n'
            f.write(fpkm_write)
    return fpkm_path

###########################################

def export_multi_gene_list(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    geneset_id = data.split(",")
    multi_geneset_path = dir_path + "/multi_geneset_list"
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    f = open(multi_geneset_path, "wb")
    for n, gi in enumerate(geneset_id):
        my_result = main_collection.find_one({'_id': ObjectId(gi)})
        if not my_result:
            bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(gi)))
        f.write(my_result["name"] + "\t")
        results = collection.find_one({"geneset_id": ObjectId(gi)})
        # id_list = []
        # for result in results:
        #     gene_id = result['gene_name']
        #     id_list.append(gene_id)
        f.write(",".join(results["gene_list"]) + "\n")
    return multi_geneset_path

def export_class_code_for_enrich(data, option_name, dir_path, bind_obj=None):
    tmp = data.split("\t")
    task_id = tmp[0]
    class_code_type = tmp[1]
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    collection = db['sg_express_class_code']
    result = collection.find_one({"task_id":task_id})
    id = result["_id"]
    class_code = os.path.join(dir_path, "%s_class_code" % task_id)
    bind_obj.logger.debug("正在导出class_code信息:%s" %(class_code))
    class_code_detail = db['sg_express_class_code_detail']
    class_code_col = db['sg_express_class_code']
    class_code_main_info = class_code_col.find_one({"_id": id})
    task_id = class_code_main_info["task_id"]
    sg_task = db["sg_task"]
    sg_task_info = sg_task.find_one({"task_id":task_id})
    if sg_task_info["is_demo"] == 2:
        demo_id = sg_task_info["demo_id"]  # sanger_21455
        class_code_demo_info = class_code_col.find_one({"task_id":demo_id})
        bind_obj.logger.info(class_code_demo_info)
        data = str(class_code_demo_info["_id"])
    bind_obj.logger.info(data)
    class_code_info = class_code_detail.find({"class_code_id":ObjectId(data),"type":class_code_type})
    with open(class_code,'w+') as f:
        header = ['seq_id','gene_name',"class_code"]
        f.write("\t".join(header)+"\n")
        for d in class_code_info:
            if type == 'gene':
                _write = d['assembly_gene_id']+"\t"+d['gene_name']+"\t" + d["class_code"] + "\n"
            if type == 'transcript':
                _write = d['assembly_trans_id']+"\t"+d['gene_name']+"\t" + d["class_code"] + "\t" + d["assembly_gene_id"] + "\n"
            f.write(_write)
    return class_code

# for venn_exp
def export_group(data, option_name, dir_path, bind_obj=None):
    from collections import OrderedDict
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    group_out = os.path.join(dir_path, option_name)
    with open(group_out, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in group_dict:
            for each in group_dict[key]:
                f.write('{}\t{}\n'.format(each, key))
    return group_out

# added  by gdq for wgcna
def sample_ids2sample_names(sample_id_list):
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    specimen_info = db['sg_specimen']
    sample_names = list()
    for each in sample_id_list:
        sample_info = specimen_info.find_one({"_id": ObjectId(each), "about_qc" : "after"})
        if sample_info:
            sample_names.append(sample_info['specimen_name'])
        else:
            bind_obj.set_error('{} not found in sg_specimen'.format(each))
    return sample_names

# added  by gdq for wgcna
def export_geneset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    from collections import OrderedDict
    import pandas as pd
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    exp_id, geneset_id = data.split(";")
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
        group_dict[each] = sample_ids2sample_names(group_dict[each])
    with open(dir_path+'/group_info.txt', 'w') as f:
        f.write('#sample\tgroup\n')
        group_id = bind_obj.sheet.option('group_id').lower()
        for g in group_dict:
            for s in group_dict[g]:
                if group_id == "all":
                    g = s
                f.write('{}\t{}\n'.format(s, g))
    samples = sample_ids2sample_names(samples)
    target_cols = OrderedDict(seq_id=1, _id=0)
    conn =  db['sg_express']
    result = conn.find_one({"_id": ObjectId(exp_id)})
    if not result:
        bind_obj.set_error("找不到制定的sg_express主表：{}".format(exp_id))
    exp_type = result['name'].split('_')[2]
    for each in samples:
        target_cols[each] = 1

    # get geneset
    if "all" not in geneset_id.lower():
        print(geneset_id)
        conn = db['sg_geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        geneset = geneset_records['gene_list']
        conn = db['sg_geneset']
        result = conn.find_one({"_id": ObjectId(geneset_id)})
        geneset_type = result['type']
    else:
        geneset = None
        geneset_type = bind_obj.sheet.option('exp_level')

    # get exp matrix
    conn = db['sg_express_detail']
    query_dict = {"express_id": ObjectId(exp_id), 'type': geneset_type, "sample_group":"sample", "value_type":exp_type}
    exp_records = conn.find(query_dict, target_cols)
    print(query_dict)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('seq_id')
    if geneset:
        exp_matrix = exp_matrix.loc[geneset,:]
    if 'refall' in geneset_id.lower():
        geneset = [x for x in exp_matrix.index if (not x.startswith('MSTRG')) and (not x.startswith('TCONS')) and (not x.startswith('XLOC'))]
        exp_matrix = exp_matrix.loc[geneset,:]
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def download_s3_file(path, to_path):
    """
    判断文件是否在对象存储上
    """
    if not to_path.startswith("/"):
        to_path = os.path.join(to_path)
    if os.path.exists(to_path):
        os.remove(to_path)
    elif os.path.exists(path):
        to_path = path
    elif exists(path):
        download(path, to_path)
    else:
        print 'file can not find'
    return to_path

# added by gdq for wgcna
def export_wgcna_exp_matrix(data, option_name, dir_path, bind_obj=None):
    import pandas as pd
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    wgcna_prepare = db['sg_wgcna_prepare']
    prepare_id = data.strip()
    prepare_main = wgcna_prepare.find_one({"_id": ObjectId(prepare_id)})
    exp_matrix = prepare_main['exp_matrix']
    if re.match(r'^\w+://\S+/.+$', exp_matrix):
        download_s3_file(exp_matrix, os.path.join(bind_obj.work_dir, "exp_matrix.txt"))
        exp_matrix = os.path.join(bind_obj.work_dir, "exp_matrix.txt")
        if not os.path.exists(os.path.join(dir_path, "exp_matrix.txt")):
            os.link(exp_matrix, os.path.join(dir_path, "exp_matrix.txt"))
    else:
        if os.path.exists(os.path.join(dir_path, "exp_matrix.txt")):
            os.remove(os.path.join(dir_path, "exp_matrix.txt"))
        os.link(exp_matrix, os.path.join(dir_path, "exp_matrix.txt"))

    target_seqs = pd.read_table(exp_matrix, header=0, index_col=0).index
    task_id = prepare_main['task_id']
    sg_task = db['sg_task']
    sg_task_info = sg_task.find_one({"task_id": task_id})
    if 'is_demo' in sg_task_info and int(sg_task_info['is_demo']) != 0:
        task_id = sg_task_info['demo_id']
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id})
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_annotation_query_detail']
    anno_type = bind_obj.sheet.option('exp_level')
    query_dict = dict(
        query_id=annot_main_id,
        anno_type=anno_type,
    )
    result_dict = dict( _id=0, gene_name=1, gene_id=1)
    if anno_type == 'transcript':
        result_dict.update({"transcript_id": 1})
    result = annot_detail.find(query_dict, result_dict)
    print(query_dict)
    print(result_dict)
    gene2name = pd.DataFrame(list(result))
    print(gene2name.head())
    if anno_type == 'transcript':
        gene2name.set_index('transcript_id', inplace=True)
        gene2name = gene2name.loc[:, ["gene_id", "gene_name"]]
    else:
        gene2name.set_index('gene_id', inplace=True)
    gene2name = gene2name.loc[list(target_seqs), :]
    gene2name.reset_index(inplace=True)
    output = os.path.join(dir_path, "seq_id2gene_name.txt")
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    gene2name = pd.read_table(output, header=0)
    gene2name.fillna(method="pad", axis=1, inplace=True)
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    return exp_matrix+';'+output

# added  by gdq for wgcna
def export_wgcna_relate_input(data, option_name, dir_path, bind_obj=None):
    import pandas as pd
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    module_id = data.strip()
    # export eigengenes
    eigengenes = db['sg_wgcna_module_eigengenes_detail']
    eigengenes_found = eigengenes.find({"module_id": ObjectId(module_id)}, {"_id": 0, "module_id":0})
    eigengenes_pd = pd.DataFrame(list(eigengenes_found))
    eigengenes_pd.set_index("module", inplace=True)
    eigengenes_path = os.path.join(dir_path, "module_eigengenes.xls")
    eigengenes_pd.to_csv(eigengenes_path, sep='\t', header=True, index=True)
    # export exp matrix
    prepare_id = db['sg_wgcna_module'].find_one({"_id": ObjectId(module_id)})['wgcna_prepare_id']
    exp_matrix = db['sg_wgcna_prepare'].find_one({"_id": ObjectId(prepare_id)})['exp_matrix']
    # export each gene's module and gene_id/gene_name info
    exp_level = bind_obj.sheet.option('exp_level')
    membership = db['sg_wgcna_module_membership_detail']
    query_dict = dict(module_id=ObjectId(module_id))
    return_dict = dict(_id=0, seq_id=1, gene_name=1, module=1, kme=1, block_id=1)
    if exp_level == 'transcript':
        return_dict.update({"gene_id": 1})
    result = membership.find(query_dict, return_dict)
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("seq_id", inplace=True)
    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    return exp_matrix + ';' + eigengenes_path + ";" + gene_annot

# added for transcription factor analysis
def get_all_pep_seq(data, option_name, dir_path, bind_obj=None):
    print  "***\n" + data
    pep_db_path, task_id = data.strip().split(',')
    if re.match(r'^\w+://\S+/.+$', pep_db_path) or re.match(r'/mnt/ilustre', pep_db_path):
        transfer = MultiFileTransfer()
        transfer.add_download(pep_db_path, bind_obj.work_dir + "/")
        transfer.perform()
        #pep_db_path = download_s3_file(pep_db_path, os.path.join(bind_obj.work_dir, "pep.db"))
        pep_db_path = os.path.join(bind_obj.work_dir, "refrna_seqs.db")
    # get transcript list
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    sg_task = db['sg_task']
    sg_task_info = sg_task.find_one({"task_id": task_id})
    if 'is_demo' in sg_task_info and int(sg_task_info['is_demo']) != 0:
        task_id = sg_task_info['demo_id']
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id})
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_annotation_query_detail']
    query_dict = dict(query_id=annot_main_id, anno_type="transcript", )
    result_dict = dict(_id=0, transcript_id=1)
    result = annot_detail.find(query_dict, result_dict)
    transcripts = set(x['transcript_id'] for x in result)

    result_path = os.path.join(dir_path, "all_pep.fa")
    import sqlite3
    conn = sqlite3.connect(pep_db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM pep")
    with open(result_path, 'w') as fw:
        # pep_id 就是转录本的id
        for pep_id, pep_seq in cursor.fetchall():
            if pep_id in transcripts:
                fw.write('>{}\n{}\n'.format(pep_id, pep_seq))
    return result_path

def get_gene_detail(data, option_name, dir_path, bind_obj=None):
    import pandas as pd
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    conn = db['sg_express_class_code_detail']
    gene_detail = conn.find(
        {"class_code_id": ObjectId(data.strip()), "is_new": {"$exists":True}},
        {"_id": 0, "gene_id": 1, "description": 1, "gene_name": 1, "transcript": 1,}
    )
    data = list()
    for each_dict in gene_detail:
        transcripts = each_dict["transcript"].split(',')
        for trans in transcripts:
            data.append(dict(
                transcript_id=trans,
                gene_id=each_dict['gene_id'],
                gene_name=each_dict['gene_name'],
                gene_desc=each_dict['description']
            ))
    result_pd = pd.DataFrame(list(data))
    result_pd.set_index("transcript_id", inplace=True)
    result_pd = result_pd.loc[:, ["gene_id", "gene_name", "gene_desc"]]
    result_pd["gene_name"][result_pd['gene_name'] == "-"] = None
    result_pd = result_pd.fillna(method="pad", axis=1)
    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    os.system(r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot

# added by gdq for tfbs predict
def export_geneid2tfid_file(data, option_name, dir_path, bind_obj=None):
    import pandas as pd
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    conn = db['sg_tf_predict_detail']
    result = conn.find({"tf_predict_id": ObjectId(data.strip())}, {"_id": 0, 'gene_id': 1, "blast_hit": 1})
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("gene_id", inplace=True)
    result_file = os.path.join(dir_path, "geneid2tfid.txt")
    result_pd.to_csv(result_file, sep='\t', header=True, index=True)
    return result_file

def export_predict_result(data, option_name, dir_path, bind_obj=None):
    import pandas as pd
    db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    conn = db['sg_tfbs_predict_detail']
    target_cols = {
        "_id":0, "tf_geneid":1, "gene_name_tf":1, "target_id":1, "gene_name_target":1,
        "corr":1, "corr_pvalue":1, "corr_padjust":1, "p-value":1, "q-value":1
    }
    result = conn.find({"tfbs_predict_id": ObjectId(data.strip()),}, target_cols,)
    result_pd = pd.DataFrame(list(result))
    result_file = os.path.join(dir_path, "tfbs_predict.xls")
    result_pd.to_csv(result_file, sep='\t', header=True, index=False)
    return result_file
