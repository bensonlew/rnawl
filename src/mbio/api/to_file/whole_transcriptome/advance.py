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
from collections import OrderedDict
import pandas as pd
from biocluster.file import getsize, exists
from biocluster.file import download
#from biocluster.api.file.lib.s3 import S3TransferManager
#from boto.s3.bucket import Bucket
import shutil
from biocluster.api.file.lib.transfer import MultiFileTransfer
import sys
import numpy as np

project_type = 'whole_transcriptome'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

'''
def download_from_s3(from_file, to_path="download/", cover=True):
    """
    s3 数据工作流， 来自框架workflow
    从s3对象存储下载数据到本地, 为了避免堵塞进程，此功能应该放置在流程最后执行。
    :param from_file: 需要下载的文件路径或文件路径, 必须是类似s3region://bucket/key写法。
    因为对象存储中没有文件夹的概念，需要下载文件夹必须使用"/"结尾，以明确表明下载的是文件夹
    :param to_path: 下载文件相对于当前工作目录的存放目录。
    当路径为"/"结尾时，表示下载文件存放在此文件夹下，否者为下载完整路径。
    当from_file为文件夹时，此参数也必须以"/"结尾。目录层级与下载的s3目录层级结构相同。
    默认情况下放置在当前模块工作目录的download目录下。
    :param cover: 对已存在的文件是否覆盖
    :return:
    """
    if re.match(r"^/|^\.\.", to_path):
        raise Exception("不能使用绝对路径或切换到其他目录!")
    if os.path.basename(to_path) == ".":
        raise Exception("目标文件不能叫\".\"!")
    target_dir = False
    if re.match(r"/$", to_path):
        target_dir = True
    work_dir = os.getcwd()
    s3transfer = S3TransferManager()
    s3transfer.base_path = work_dir
    s3transfer.overwrite = cover
    m = re.match(r"^([\w\-]+)://([\w\-]+)/(.*)$", from_file)
    if not m:
        raise Exception("下载路径%s格式不正确!" % from_file)
    else:
        region = m.group(1)
        bucket_name = m.group(2)
        key_name = m.group(3)
        if re.match(r"/$", key_name):
            if not target_dir:
                raise Exception("下载文件为文件夹时，源路径%s也必须为文件夹,以\"/\"结尾!" % to_path)
            conn = s3transfer.config.get_rgw_conn(region, bucket_name)
            bucket = Bucket(connection=conn, name=bucket_name)
            for key in bucket.list(prefix=key_name):
                source = os.path.join(from_file, key.name)
                target = os.path.join(target_dir, os.path.relpath(key.name, key_name))
                s3transfer.add(source, target)
        else:
            if not target_dir:  # 处理已存在文件的情况
                target = os.path.join(work_dir, to_path)
                if os.path.exists(target):
                    if cover:
                        if os.path.isdir(target):
                            shutil.rmtree(target)
                        else:
                            os.remove(target)
                    else:
                        raise Exception("目标文件夹%s已经存在!" % target)
            else:
                target = os.path.join(work_dir, to_path, os.path.basename(key_name))
            s3transfer.add(from_file, target)
    s3transfer.wait()
'''

def export_gene_list(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "%s_gene.list" % option_name)
    bind_obj.logger.debug("正在导出基因集")
    collection = db['geneset_detail']
    main_collection = db['geneset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:%s在geneset中未找到！", variables=(ObjectId(data)), code="51008132")
    results = collection.find_one({"geneset_id": ObjectId(data)})
    with open(gene_list_path, "wb") as f:
        gene_list = results["seq_list"]
        for gene_id in gene_list:
            f.write(gene_id + "\n")
    return gene_list_path


def export_circ_detail(data, option_name, dir_path, bind_obj=None):
    bind_obj.logger.debug("circ详情")
    collection = db['circrna_identify_detail']
    main_collection = db['circrna_identify']
    my_result = main_collection.find_one({'task_id': data})
    main_id = my_result["main_id"]
    if not my_result:
        return None
        # bind_obj.set_error("意外错误，未找到circ rna相关数据！")
    results = collection.find({"detail_id": main_id})

    circ_detail = dir_path + '/circ_detail'
    with open(circ_detail, "wb") as w:
        w.write("circrna_id\thost_gene_id\n")
        for result in results:
            w.write("{}\t{}\n".format(result['circrna_id'], result['host_gene_id']))
    return circ_detail

def export_mrna_detail(data, option_name, dir_path, bind_obj=None):
    task_id = data
    exp_collection = db['exp']
    exp_detail_collection = db['exp_detail']
    main_id =  exp_collection.find_one({"task_id": task_id, "level": "T"})["main_id"]
    genes_detail = exp_detail_collection.find({"exp_id": main_id, 'category': "mRNA", 'level': "T"})
    with open(dir_path + "/mrna_annot", "w") as w:
        w.write("gene_id\ttranscript_id\tis_gene\tgene_name\tlength\tdescription\n")

        for gene_dict in genes_detail:
            w.write("\t".join([
                gene_dict.get("gene_id", ""),
                gene_dict.get("transcript_id", ""),
                "yes",
                gene_dict.get("gene_name", ""),
                str(gene_dict.get("length", 0)),
                gene_dict.get("description", "")
            ]) + "\n")
    return dir_path + "/mrna_annot"

def get_mirna_seq(task_id, category, kind, set_list, dir_path, bind_obj=None):
    genes_collection = db['genes']
    genes_detail_collection = db['genes_detail']
    main_id =  genes_collection.find_one({"task_id": task_id})["main_id"]
    if set_list == "all":
        seq_select_list = genes_detail_collection.find({'genes_id':ObjectId(main_id), 'category': category, 'level': "T"})
    else:
        seq_select_list = genes_detail_collection.find({'genes_id':ObjectId(main_id), 'category': category, 'level': "T", 'transcript_id': {"$in": set_list}})

    seq_path = dir_path + "/" +  category + "_" + kind + ".fa"
    gene_set = list()
    with open(seq_path, "wb") as w:
        for seq in seq_select_list:
            if seq['transcript_id'] in gene_set:
                pass
            else:
                w.write(">{}\n{}\n".format(seq['transcript_id'], seq['mirna_seq']))
                gene_set.append(seq['transcript_id'])
    return seq_path

def get_longrna_seq(task_id, category, kind, set_list, dir_path, bind_obj=None):
    task_collection = db['task']
    try:
        out_put = task_collection.find_one({"task_id": task_id})["output"]
        if category == "mRNA":
            fasta_dir = out_put + "/mrna/02_Basic_Analysis/03_Assemble/Sequence/all_transcript.fa"
        else:

            fasta_dir = out_put + "/mrna/02_Basic_Analysis/03_Assemble/Sequence/all_transcript.fa"
            fasta_dir_lnc = out_put + "/lncrna/01_lncRNA_Analysis/01_Known_lncRNA/known_lncrna.fa"
        # long_dir = task_collection.find_one({"task_id": task_id})["sub_output"]["long"]
    except:
        bind_obj.set_error("无法找到序列文件目录")

    # fasta_dir = long_dir + "/large_gush/merge_known_new/known_and_new.fa"
    download_tmp = dir_path + "/"+ '{}tmp.fasta'.format(np.random.randint(10000))
    download_s3_file(fasta_dir, download_tmp)
    if category == "lncRNA":
        download_tmp2 = dir_path + "/"+ '{}tmp2.fasta'.format(np.random.randint(10000))
        download_s3_file(fasta_dir_lnc, download_tmp2)

    from Bio import SeqIO
    seq_path = dir_path + "/" +  category + "_" + kind + ".fa"
    gene_set = set(set_list)
    with open(seq_path, "wb") as w:
        for seq in SeqIO.parse(download_tmp, "fasta"):
            if set_list == "all" or seq.id in gene_set:
                w.write(">{}\n{}\n".format(seq.id, seq.seq))
        if category == "lncRNA":
            for seq in SeqIO.parse(download_tmp2, "fasta"):
                if set_list == "all" or seq.id in gene_set:
                    w.write(">{}\n{}\n".format(seq.id, seq.seq))
    return seq_path


def get_circrna_seq(task_id, category, kind, set_list, dir_path, bind_obj=None):
    task_collection = db['task']
    try:
        out_put = task_collection.find_one({"task_id": task_id})["output"]

        fasta_dir = out_put + "/circrna/03_circRNA_Predict/circRNA.fa"
        '''
        # sub_output = task_collection.find_one({"task_id": task_id})["sub_output"]
        fasta_dir = sub_output["long"] + "/circ_brush/circrna.fasta"
        if "circle" in sub_output:
            fasta_dir = sub_output["circle"] + "/circ_brush/circrna.fasta"
        else:
            fasta_dir = sub_output["long"] + "/circ_brush/circrna.fasta"
        '''
    except:
        bind_obj.set_error("无法找到序列文件目录")



    download_tmp = dir_path + "/" + '{}tmp.fasta'.format(np.random.randint(10000))
    download_s3_file(fasta_dir, download_tmp)
    from Bio import SeqIO
    seq_path = dir_path + "/" +  category + "_" + kind + ".fa"
    gene_set = set(set_list)
    if set_list == "all":
        return download_tmp
    else:
        with open(seq_path, "wb") as w:
            for seq in SeqIO.parse(download_tmp, "fasta"):
                if set_list == "all" or seq.id in gene_set:
                    w.write(">{}\n{}\n".format(seq.id, seq.seq))
        return seq_path


def get_seq(data, option_name, dir_path, bind_obj=None):
    task_id, category, kind, gene_set = data.split(";")
    print "data is {}".format(data)
    if gene_set.lower() != "all":
        collection = db['geneset_detail']
        main_collection = db['geneset']
        my_result = collection.find_one({'geneset_id': ObjectId(gene_set)})

        gene_set_list = zip(my_result['seq_list'], my_result['kind_list'], my_result['category_list'])

        if kind == "all":
            set_list = [gene[0] for gene in gene_set_list if gene[2] == category]
        else:
            set_list = [gene[0] for gene in gene_set_list if gene[2] == category and gene[1] == kind]
    else:
        # 基因集为all时 exp_detail 中判断基因类型
        main_collection = db['exp']
        my_result = main_collection.find_one({'task_id': task_id, 'level': 'T'})

        collection = db['exp_detail']
        results = collection.find({"exp_id": ObjectId(my_result["main_id"]), "category": category, "level": "T", "kind": kind})
        set_list = [result["transcript_id"] for result in results]

    if category == "miRNA":
        seq_path = get_mirna_seq(task_id, category, kind, set_list, dir_path, bind_obj=bind_obj)
    elif category == "mRNA":
        seq_path = get_longrna_seq(task_id, category, kind, set_list, dir_path, bind_obj=bind_obj)
    elif category == "lncRNA":
        seq_path = get_longrna_seq(task_id, category, kind, set_list, dir_path, bind_obj=bind_obj)
    elif category == "circRNA":
        seq_path = get_circrna_seq(task_id, category, kind, set_list, dir_path, bind_obj=bind_obj)
    print "get seq_path {}".format(seq_path)
    return seq_path


def export_go_list(data, option_name, dir_path, bind_obj=None):
    go_list_path = os.path.join(dir_path, "GO.list")
    bind_obj.logger.debug("正在导出%sgo列表:%s" % (option_name, go_list_path))
    geneset_collection = db["geneset"]
    # task_id = geneset_collection.find_one({"main_id": ObjectId(data)})["task_id"]
    my_result = db["annotation_go"].find_one({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    geneset_type= bind_obj.sheet.option('geneset_type')
    go_id = my_result["main_id"]
    if not my_result:
        bind_obj.set_error("意外错误，annotation_go_id:%s在annotation_go中未找到！", variables=(go_id), code="51008133")
    collection = db["annotation_go_list"]
    results = collection.find({"go_id": ObjectId(go_id)})
    one_record = collection.find_one({"go_id": ObjectId(go_id)})
    if not one_record:
        bind_obj.set_error("生成gos_list出错：annotation_id:%s在annotation_gos_list中未找到！", variables=(ObjectId(go_id)), code="51008134")
    with open(go_list_path, "wb") as w:
        for result in results:
            gene_id = result["gene_id"]
            go_list = result["gos_list"]
            go_anno_type = result["anno_type"]
            if go_anno_type == geneset_type:
                w.write(gene_id + "\t" + go_list + "\n")
    return go_list_path

# 通过基因集id到geneset去获取task_id，type对应的信息，然后到annotation_kegg去查找注释信息，然后导表；页面上的选择限制2个geneset_id的类型
# 肯定是一样，所以任意选择一个geneset_id过来就可以，所以接口那里随便选择了一个
# 这里的one_record是多余的，直接find的结果就可以判断，既然写了就写了吧
def export_kegg_table(data, option_name, dir_path, bind_obj=None):
    kegg_path = os.path.join(dir_path, 'gene_kegg_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg_path))
    geneset_collection = db["geneset"]
    bind_obj.logger.debug(data)
    geneset_result = geneset_collection.find_one({"main_id": ObjectId(data)})
    task_id = geneset_result["task_id"]
    bind_obj.logger.debug("ttttttt")
    bind_obj.logger.debug(task_id)
    geneset_type = geneset_result["level"]
    my_result = db["annotation_kegg"].find({"task_id": bind_obj.sheet.option("task_id")})
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\n')
        for main_table in my_result:
            kegg_id = main_table["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not my_result:
                bind_obj.set_error("意外错误，annotation_kegg_id:%s在annotation_kegg中未找到！", variables=(kegg_id), code="51008135")
            bind_obj.logger.debug("{} 搜索条件 {} ".format("annotation_kegg_table", {'kegg_id': kegg_id, 'anno_type': geneset_type}))
            results = db['annotation_kegg_table'].find({'kegg_id': kegg_id, 'anno_type': geneset_type})
            one_record = db['annotation_kegg_table'].find_one({'kegg_id': kegg_id, 'anno_type': geneset_type})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:%s在annotation_kegg_table中未找到！", variables=(ObjectId(kegg_id)), code="51008136")
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['ko_id'] + result['transcript_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths']))
    return kegg_path

def export_kegg_level_table(data, option_name, dir_path, bind_obj=None):
    kegg__level_path = os.path.join(dir_path, 'gene_kegg_level_table.xls')
    bind_obj.logger.debug("正在导出参数%s的kegg_table文件，路径:%s" % (option_name, kegg__level_path))
    geneset_collection = db["geneset"]
    bind_obj.logger.debug(data)
    geneset_result = geneset_collection.find_one({"main_id": ObjectId(data)})
    task_id = geneset_result["task_id"]
    bind_obj.logger.debug(task_id)
    my_result = db["annotation_kegg"].find({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option("type")})
    with open(kegg__level_path, 'wb') as w:
        w.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
        for i in my_result:
            kegg_id = i["main_id"]
            bind_obj.logger.debug(kegg_id)
            if not kegg_id:
                bind_obj.set_error("意外错误，annotation_kegg_id:%s在annotation_kegg中未找到！", variables=(kegg_id), code="51008137")
            results = db["annotation_kegg_level"].find({"kegg_id": kegg_id, "seq_type": "all", "anno_type": bind_obj.sheet.option("geneset_type")})
            one_record = db['annotation_kegg_level'].find_one({'kegg_id': kegg_id, "seq_type": "all", 'anno_type': bind_obj.sheet.option("geneset_type")})
            if not one_record:
                results = db["annotation_kegg_level"].find(
                    {"kegg_id": kegg_id, "seq_type": "ref", "anno_type": bind_obj.sheet.option("geneset_type")})
                one_record = db['annotation_kegg_level'].find_one(
                    {'kegg_id': kegg_id, "seq_type": "ref", 'anno_type': bind_obj.sheet.option("geneset_type")})
            if not one_record:
                bind_obj.set_error("生成kegg_table出错：kegg_id:%s在annotation_kegg_level中未找到！", variables=(ObjectId(kegg_id)), code="51008138")
            for result in results:
                if 'hyperlink' not in result:
                    bind_obj.logger.debug(result['pathway_id'] + result['graph_id'] + '-> no hyperlink')
                    result['hyperlink'] = 'None'
                w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '', result['number_of_seqs'], result['pathway_definition'],
                result['first_category'], result['anno_type'], result['hyperlink'], result['seq_list'], '', result['second_category']))
    return kegg__level_path

def export_kegg_pdf(data, option_name, dir_path, bind_obj=None):
    fs = gridfs.GridFS(db)
    annotation_collection = db["annotation_kegg"]
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    main_id = annotation_collection.find_one({"task_id":task_id, "type": bind_obj.sheet.option("type")})["main_id"]
    kegg_level_collection = db["annotation_kegg_level"]
    results = kegg_level_collection.find({"kegg_id":main_id, "anno_type":anno_type})
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
    all_list = os.path.join(dir_path, "all_gene.list")
    bind_obj.logger.debug("正在导出所有背景基因{}".format(all_list))
    collection = db['exp_detail']
    exp_collection = db['exp']
    main_collection = db['geneset']
    my_result = main_collection.find_one({"main_id": ObjectId(data)})
    task_id = my_result["task_id"]
    geneset_type = my_result["level"]
    bind_obj.logger.debug(task_id)
    exp_result = exp_collection.find_one({'task_id': bind_obj.sheet.option('task_id'), 'exp_level': geneset_type})
    if not exp_result:
        bind_obj.set_error("意外错误，task_id:%s的背景基因在geneset中未找到！", variables=(data), code="51008139")
    exp_id = exp_result["main_id"]
    results = collection.find({"exp_id": ObjectId(exp_id)})
    with open(all_list, "wb") as f:
        for result in results:
            gene_id = result['seq_id']
            f.write(gene_id + "\n")
    return all_list


def get_gene_link(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip()
    annot_table = db['genes']
    annot_main = annot_table.find_one({"task_id": task_id})
    if not annot_main:
        return "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/whole_rna/test.trait.txt.checked.txt"
        # Exception("Not Found in annotation_query by query {}".format(task_id))
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    lnc_detail = db['genes_detail']
    query_dict = dict(main_id=annot_main_id)
    result_dict = dict(_id=0, gene_id=1, gene_ensembl=1)
    result = lnc_detail.find(query_dict, result_dict)
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("gene_id", inplace=True)
    result_pd = result_pd.loc[:, ["gene_id", "gene_ensembl"]]
    # result_pd.columns = ["gene_id", "gene_lnc"]
    gene_lnc = os.path.join(dir_path, "gene_lnc.xls")
    result_pd.to_csv(gene_lnc, sep='\t', header=True, index=True)
    return gene_lnc


def export_cog_class(data, option_name, dir_path, bind_obj=None):
    cog_path = os.path.join(dir_path, 'cog_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(cog_path))
    genesets, table_title, task_id, geneset_type = get_geneset_detail(data, bind_obj)
    cog_collection = db["annotation_cog"]
    cog_detail_collection = db["annotation_cog_detail"]
    cog_id = cog_collection.find_one({"task_id": bind_obj.sheet.option('task_id'), 'type': bind_obj.sheet.option('type')})["main_id"]
    print("cog_id:", cog_id, geneset_type,  )
    cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type})
    # 更改mongo查询是否为空的判断
    num = 0
    for i in cog_results:
        num += 1
    if num == 0:
        cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "ref", 'anno_type': geneset_type})
    else:
        # cog_results循环完了需要重新赋值，否则为空
        cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type})
    # if cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type}):
    #     cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "all", 'anno_type': geneset_type})
    # else:
    #     cog_results = cog_detail_collection.find({'cog_id': cog_id, "seq_type": "ref", 'anno_type': geneset_type})
    new_table_title = []
    for tt in table_title:
        # new_tt = [tt + "_COG", tt + "_NOG", tt + "_KOG", tt + "_COG_list", tt + "_NOG_list", tt + "_KOG_list"]
        new_tt = [tt + "_COG", tt + "_COG_list"]
        new_table_title = new_table_title + new_tt
    bind_obj.logger.debug(table_title)
    with open(cog_path, "wb") as w:
        w.write("Type\tFunctional Categoris\t" + "\t".join(new_table_title) + "\n")
        for cr in cog_results:
            # kog_list = set([])
            # nog_list = set(cr["nog_list"].split(";") if cr["nog_list"] else [])
            cog_list = set(cr["cog_list"].split(";") if cr["cog_list"] else [])
            write_line = {}
            for gt in genesets:
                # kog_count = list(kog_list & genesets[gt][1])
                # nog_count = list(nog_list & genesets[gt][1])
                print gt
                print genesets[gt][1]
                print cog_list
                cog_count = list(cog_list & genesets[gt][1])
                # if not len(kog_count) + len(nog_count) + len(cog_count) == 0:
                #     write_line[gt] = [str(len(cog_count)), str(len(nog_count)), str(len(kog_count)), ";".join(cog_count), ";".join(nog_count), ";".join(kog_count)]
                if not len(cog_count) == 0:
                    write_line[gt] = [str(len(cog_count)), ";".join(cog_count)]
            print('here:',  genesets[gt][1])

            if len(write_line) > 0:
                w.write("{}\t{}\t".format(cr["type"], cr["function_categories"]))
                for tt in table_title:
                    w.write("\t".join(write_line[tt]) + "\t") if tt in write_line else w.write("0\tnone\t")
                w.write("\n")
    return cog_path


def get_geneset_detail(data, bind_obj):
    geneset_collection = db["geneset"]
    genesets = {}
    names = []
    task_id = ""
    geneset_type = bind_obj.sheet.option('geneset_type')
    for geneset_id in data.split(","):
        geneset_result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not geneset_result:
            bind_obj.set_error("意外错误:未找到基因集_id为%s的基因集信息", variables=(geneset_id), code="51008140")
        task_id = geneset_result["task_id"]
        geneset_type = geneset_result["level"]
        geneset_name = geneset_result["name"]
        names.append(geneset_name)
        genesets[geneset_name] = [geneset_type]
        collection = db['geneset_detail']
        results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
        geneset_names = set(results["seq_list"])
        genesets[geneset_name].append(geneset_names)
    #print genesets
    return genesets, names, task_id, geneset_type


def export_go_class(data, option_name, dir_path, bind_obj=None):
    go_path = os.path.join(dir_path, 'go_class_table.xls')
    bind_obj.logger.debug("正在导出{}".format(go_path))
    genesets, names, task_id, seq_type = get_geneset_detail(data, bind_obj)
    go_collection = db["annotation_go"]
    go_level_collection = db["annotation_go_detail"]
    go_id = go_collection.find_one({"task_id": bind_obj.sheet.option('task_id'), "type": bind_obj.sheet.option('type')})["main_id"]
    one_record = go_level_collection.find_one({'go_id': go_id, "level": 2, "anno_type": seq_type})
    if not one_record:
        bind_obj.set_error("意外错误:未找到go_id为%s的基因集信息", variables=(go_id), code="51008141")
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
            if go_level_collection.find_one({'go_id': go_id, "seq_type": "all", "level": 2, "anno_type": seq_type}) :
                go_results = go_level_collection.find({'go_id': go_id, "seq_type": "all", "level": 2, "anno_type": seq_type})
            else:
                go_results = go_level_collection.find(
                    {'go_id': go_id, "seq_type": "ref", "level": 2, "anno_type": seq_type})
            for gr in go_results:
                if gr["goterm"] == item:
                    seq_list = set(gr["seq_list"].split(";"))
                    write_line = {}
                    for gt in genesets:
                        total_gene_num = len(genesets[gt][1])
                        go_count = list(seq_list & genesets[gt][1])
                        if not len(go_count) == 0:
                            write_line[gt] = str(len(go_count)) + "\t" + str(len(go_count)/total_gene_num) + \
                                             "(" + str(len(go_count)) + "/" + str(total_gene_num) + ")" + "\t" + ";".join(go_count)
                    if len(write_line):
                        w.write("{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"]))
                        for tt in genesets:
                            w.write(write_line[tt] + "\t") if tt in write_line else w.write("0\t0\tnone\t")
                        w.write("\n")
    return go_path


# ############表达量部分

def export_group_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    按分组的详细信息获取group表
    使用时确保你的workflow的option里group_detal这个字段
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的GROUP表格为文件，路径:%s" % (option_name, file_path))
    if data in ["all", "All", "ALL"]:
        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
        return file_path
    data = _get_objectid(data)
    group_detail = bind_obj.sheet.option('group_detail')  #另传字段
    group_table = db['specimen_group']
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            bind_obj.set_error("生成group表失败，传入的%s不是一个字典或者是字典对应的字符串", variables=(option_name), code="51008142")
    if not isinstance(table_dict, dict):
        bind_obj.set_error("生成group表失败，传入的%s不是一个字典或者是字典对应的字符串", variables=(option_name), code="51008143")
    group_schema = group_table.find_one({"main_id": ObjectId(data)})
    if not group_schema:
        bind_obj.set_error("无法根据传入的group_id:%s在specimen_group_compare表里找到相应的记录", variables=(data), code="51008144")
    schema_name = re.sub("\s", "_", group_schema["group_name"])  # 将分组方案名的空格替换成下划线
    with open(file_path, "wb") as f:
        f.write("#sample\t" + "group" + "\n")
    sample_table_name = 'specimen'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"main_id": ObjectId(sp_id)})
                if not sp:
                    bind_obj.set_error("group_detal中的样本_id:%s在样本表%s中未找到", variables=(sp_id, sample_table_name), code="51008145")
                else:
                    sp_name = sp["specimen_name"]
                f.write("{}\t{}\n".format(sp_name, k))
    return file_path


def _get_objectid(data):
    if not isinstance(data, ObjectId):
        if not isinstance(data, StringTypes):
            pass
            # bind_obj.set_error("%s不为ObjectId类型或者其对应的字符串", variables=(data), code="51008146")
        else:
            try:
                data = ObjectId(data)
            except:
                pass
                # bind_obj.set_error("%s不为ObjectId类型或者其对应的字符串", variables=(data), code="51008147")
    return data


def export_control_file(data, option_name, dir_path, bind_obj=None):  #此函数待定 不一定对
    file_path = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出计数矩阵:%s" % file_path)
    collection = db['specimen_group_compare']
    result = collection.find_one({'main_id': ObjectId(data)})
    if not result:
        bind_obj.set_error("意外错误，control_id:%s在specimen_group_compare中未找到！", variables=(ObjectId(data)), code="51008148")
    group_id = result['specimen_group_id']
    if group_id not in ['all', 'All', 'ALL']:
        """检查group_id的信息"""
        if isinstance(group_id, types.StringTypes):
            group_id = ObjectId(group_id)
        group_coll = db['specimen_group']
        g_result = group_coll.find_one({'main_id': group_id})
        if not g_result:
            bind_obj.set_error("意外错误，control_file的group_id:%s在specimen_group中未找到！", variables=(group_id), code="51008149")
    control_detail = json.loads(result['compare_names'])
    with open(file_path, 'wb') as w:
        w.write('#control\t{}\n'.format('group'))
        for i in control_detail:    #此处需要修改, 可能会有错误
            # w.write('{}\t{}\n'.format(i.keys()[0], i.values()[0]))
            control_other = i.split("|")
            w.write('{}\t{}\n'.format(control_other[0], control_other[1]))
    return file_path


def _get_gene_id(geneset, geneset_detail, _id):
    try:
        results = geneset_detail.find_one({"geneset_id":ObjectId(_id)})
        seq_id = results['seq_list']
    except Exception:
        pass
        # bind_obj.set_error("{}在geneset_detail表中没有找到!", code="51008150")
    try:
        my_result = geneset.find_one({"main_id":ObjectId(_id)})
        _name = my_result['name']
    except Exception:
        pass
        # bind_obj.set_error("{}在geneset表中没有找到!", code="51008151")
    return seq_id, _name

# 根据"add_info":geneset_info['task_id'] + "\t" + data.geneset_type，就是geneset主表的task_id和页面传过来的是annotation还是annotation1
# 也就是是origin还是latest
def export_add_info(data, option_name,dir_path, bind_obj=None):
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出add_info信息")
    col = db["annotation_kegg"]
    result = col.find_one({"task_id":task_id, "type": bind_obj.sheet.option('type')})
    insert_id = result["main_id"]
    col = db["annotation_kegg_level"]
    results = col.find({"kegg_id":insert_id, "anno_type":anno_type})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info


def export_multi_gene_list(data, option_name, dir_path, bind_obj=None):
    geneset_id = data.split(",")
    multi_geneset_path = dir_path + "/multi_geneset_list"
    collection = db['geneset_detail']
    main_collection = db['geneset']
    f = open(multi_geneset_path, "wb")
    for n, gi in enumerate(geneset_id):
        my_result = main_collection.find_one({'main_id': ObjectId(gi)})
        if not my_result:
            bind_obj.set_error("意外错误，geneset_id:%s在geneset中未找到！", variables=(ObjectId(gi)), code="51008152")
        f.write(my_result["name"] + "\t")
        results = collection.find_one({"geneset_id": ObjectId(gi)})
        f.write(",".join(results["seq_list"]) + "\n")
    return multi_geneset_path


# ---------------------表达量相关gdq----------------------------------

def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} - {}'.format(k, v))
    else:
        pass

def export_exp_matrix(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    conn = db['exp_detail']
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    exp_records = conn.find({"exp_id": ObjectId(data)}, target_cols)
    if 'type' in bind_obj.sheet.options():
        print("parameter 'type' found, and we will decide ref or new ")
        new_exp_records = list()
        if bind_obj.sheet.option('type') == "ref":
            for record in exp_records:
                if not record['seq_id'].startswith(('MSTRG', 'TCONS', 'XLOC')):
                    new_exp_records.append(record)
            exp_records = new_exp_records
        elif bind_obj.sheet.option('type') == 'new':
            for record in exp_records:
                if record['seq_id'].startswith(('MSTRG', 'TCONS', 'XLOC')):
                    new_exp_records.append(record)
            exp_records = new_exp_records
        else:
            pass
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('seq_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output


def export_group(data, option_name, dir_path, bind_obj=None):
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    group_out = os.path.join(dir_path, option_name)
    with open(group_out, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in group_dict:
            for each in group_dict[key]:
                f.write('{}\t{}\n'.format(each, key))
    return group_out


def export_compare(data, option_name, dir_path, bind_obj=None):
    conn = db['specimen_group_compare']
    result = conn.find_one({'main_id': ObjectId(data)})
    if not result:
        bind_obj.set_error("control_id: %s is not found in specimen_group_compare", variables=(ObjectId(data)), code="51008153")
    cmp_info = json.loads(result['compare_names'])
    cmp_out = os.path.join(dir_path, option_name)
    with open(cmp_out, 'w') as f:
        f.write('#ctrl\ttest\n')
        for each in cmp_info:
            f.write(each.replace('|', '\t') + '\n')
    return cmp_out
# ---------------------表达量相关----------------------------------


# ---------------------基因集及相关gdq----------------------------------
def export_geneset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    print "data is {}".format(data)
    exp_id, geneset_id, category, level = data.split(";")
    output = os.path.join(dir_path, option_name)

    category_list = category.split(",")
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    # export group info
    with open(dir_path+'/group_info.txt', 'w') as f:
        f.write('#sample\tgroup\n')
        group_id = bind_obj.sheet.option('group_id').lower()
        for g in group_dict:
            for s in group_dict[g]:
                if group_id == "all":
                    g = s
                f.write('{}\t{}\n'.format(s, g))
    # decide output columns
    if level == "G":
        target_cols = OrderedDict(gene_id=1, _id=0)
    else:
        target_cols = OrderedDict(transcript_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    # get geneset
    if "all" not in geneset_id.lower():
        conn = db['geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        if not geneset_records:
            bind_obj.set_error('geneset not found by query: %s', variables=(geneset_id), code="51008154")
        geneset = geneset_records['seq_list']
    else:
        geneset = list()

    # get all exp matrix
    conn = db['exp_detail']

    if "refall" in geneset_id.lower():
        exp_records = conn.find({"exp_id": ObjectId(exp_id), "category": {"$in": category_list}, 'kind': "ref"}, target_cols)
        bind_obj.logger.debug('exp_id {} category {} ref_all'.format(ObjectId(exp_id), category_list))
    else:
        exp_records = conn.find({"exp_id": ObjectId(exp_id), "category": {"$in": category_list}}, target_cols)
        bind_obj.logger.debug('exp_id {} category {} all'.format(ObjectId(exp_id), category_list))

    exp_matrix = pd.DataFrame(list(exp_records))
    '''
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: %s', variables=(exp_id), code="51008155")
    '''
    exp_matrix.to_csv(output + "tmp", sep='\t', header=True, index=True)
    if exp_matrix.shape[0] == 0:
        pass
    else:
        if level == "G":
            exp_matrix = exp_matrix.set_index('gene_id')
        else:
            exp_matrix = exp_matrix.set_index('transcript_id')
        bind_obj.logger.debug('gene set {}'.format(geneset))
        if len(geneset) != 0:
            geneset_choose = set(exp_matrix.index).intersection(set(geneset))
            exp_matrix = exp_matrix.loc[geneset_choose, :]


    if bind_obj.sheet.name.endswith("gsea"):
        exp_matrix.index.name = 'seq_id'

    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_geneset_mitarget(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    target_id, geneset_id, targetset_id = data.split(";")


    if "all" not in geneset_id.lower():
        conn = db['geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        if not geneset_records:
            bind_obj.set_error('geneset not found by query: %s', variables=(geneset_id), code="51008154")
        geneset = geneset_records['seq_list']
    else:
        geneset = list()


    if "all" not in targetset_id.lower():
        conn = db['geneset_detail']
        targetset_records = conn.find_one({"geneset_id": ObjectId(targetset_id)})
        if not targetset_records:
            bind_obj.set_error('targetset not found by query: %s', variables=(targetset_id), code="51008154")
        targetset = targetset_records['seq_list']
    else:
        targetset = list()

    conn = db['target_mi_detail']
    target_cols = OrderedDict(small_rna=1, target=1, gene=1, name=1, category=1, mirtarbase=1, _id=0)
    search_dict = {"target_id": ObjectId(target_id)}

    if "all" not in geneset_id.lower():
        search_dict.update({"small_rna": {"$in": geneset}})
    if "all" not in targetset_id.lower():
        search_dict.update({"target": {"$in": targetset}})

    retain_type = bind_obj.sheet.option("retain_type")
    if retain_type == "lnc2m":
        search_dict.update({"category": {"$in": ["lncRNA", "mRNA"]}})
    if retain_type == "circ2m":
        search_dict.update({"category": {"$in": ["circRNA", "mRNA"]}})

    print "search dict is {}".format(search_dict)
    target_records = conn.find(search_dict, target_cols)

    if target_records.count() == 0:
        bind_obj.set_error("未找到相关基因集的靶基因数据，请重新运行")
    output = dir_path + "/target.xls"
    col_list = [key for key in target_cols.keys() if target_cols[key] == 1]
    with open(output, 'w') as fo:
        fo.write("\t".join(col_list) + "\n")
        for re in target_records:
            fo.write("\t".join([re.get(key, "") for key in col_list]) + "\n")
    '''

    target_matrix = pd.DataFrame(list(target_records))

    output = dir_path + "/target.xls"
    target_matrix.to_csv(output, sep='\t', header=True, index=False)
    '''
    print('success to export target matrix')
    return output



def export_kegg_enrich_info_ids(data,option_name,dir_path,bind_obj=None):
    collection=db["geneset_kegg_enrich_detail"]
    if 'pathway_ids' in bind_obj.sheet.options():
        print("parameter 'pathway_ids' found, and we will choose ids for id_list ")
    ids=bind_obj.sheet.option('pathway_ids').split(",")
    target_cols = OrderedDict(id=1, _id=0,first_category=1,second_category=1,seq_list=1)
    enrich_records = collection.find({"kegg_enrich_id": ObjectId(data)}, target_cols)
    new_enrich_records = list()
    for record in enrich_records:
        if record['id'] in ids :
            new_enrich_records.append(record)
    enrich_records = new_enrich_records
    if len(enrich_records) == 0:
        bind_obj.set_error("filted_records is empty", code="51008156")
    enrich_matrix = pd.DataFrame(list(enrich_records))
    enrich_matrix = enrich_matrix.set_index('id')
    output = os.path.join(dir_path, option_name)
    enrich_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export geneset kegg enrich detail matrix')
    return output

def export_kegg_enrich_info_filter(data,option_name,dir_path,bind_obj=None):
    collection = db["geneset_kegg_enrich_detail"]
    target_cols = OrderedDict(id=1, _id=0, first_category=1, second_category=1, seq_list=1)
    new_enrich_records = list()
    if bind_obj.sheet.option("stat_level"):
        print("parameter 'stat_level' found, and we will filter data by {}").format(bind_obj.sheet.option("stat_level"))
    if bind_obj.sheet.option("stat_threshold_value") and bind_obj.sheet.option("stat_numbers_value"):
        print("stat_threshold_value is {},stat_numbers_value is {}").format(bind_obj.sheet.option("stat_threshold_value"),bind_obj.sheet.option("stat_numbers_value"))
    v1=float(bind_obj.sheet.option("stat_threshold_value"))
    v2=int(bind_obj.sheet.option("stat_numbers_value"))
    if bind_obj.sheet.option("stat_level")=="pvalue":
        target_cols["pvalue"] = 1
        enrich_records = collection.find({"kegg_enrich_id": ObjectId(data)}, target_cols)
        for record in enrich_records:
           if  float(record["pvalue"]) <= v1:
               new_enrich_records.append(record)
        enrich_records=new_enrich_records
    else:
        target_cols["corrected_pvalue"] = 1
        enrich_records = collection.find({"kegg_enrich_id": ObjectId(data)}, target_cols)
        for record in enrich_records:
            if record["corrected_pvalue"] <= v1:
                new_enrich_records.append(record)
        enrich_records = new_enrich_records
    if len(enrich_records) < v2 :
        bind_obj.set_error("filted_records less than %s ", variables=(bind_obj.sheet.option("stat_numbers_value")), code="51008157")
    enrich_matrix = pd.DataFrame(list(enrich_records))
    enrich_matrix = enrich_matrix.set_index('id')
    enrich_matrix=enrich_matrix.sort_values(bind_obj.sheet.option("stat_level"),ascending=True)[:v2]
    output = os.path.join(dir_path, option_name)
    enrich_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export geneset enrich filterd matrix')
    return output




def export_geneset_from_query(data, option_name, dir_path, bind_obj=None):
    collection = db['annotation_query_detail']
    main_collection = db['annotation_query']
    my_result = main_collection.find_one({'task_id': data})
    if not my_result:
        bind_obj.set_error("意外错误，task_id:%s在annotation_query中未找到！", variables=(data), code="51008158")
    if 'main_id' in my_result:
        query_id = my_result["main_id"]
    else:
        query_id = my_result['_id']
    results = collection.find({"query_id": ObjectId(query_id)})
    output = os.path.join(dir_path, "all_genes_transcripts.txt")
    with open(output, "w") as f:
        f.write('transcript_id\tgene_id\n')
        for result in results:
            transcript_id=''
            if u'transcript_id' in result:
                transcript_id=result['transcript_id']
            gene_id=''
            if u'gene_id' in result:
                gene_id = result['gene_id']
            f.write('{}\t{}\n'.format(transcript_id, gene_id))
    return output


def export_wgcna_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数仅为wgcna module分析服务。
    该函数会返回wgcna_prepare得到的exp_matrix路径，
    同时还会通过查询annotation_query表导出geneid和genename的关系文件:seq_id2gene_name.txt
    """
    wgcna_prepare = db['wgcna_prepare']
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
    annot_table = db['exp']
    level = bind_obj.sheet.option('level')
    try:
        if level.upper() == 'T':
            annot_main = annot_table.find_one({"task_id": task_id, "level": "T"})
        else:
            annot_main = annot_table.find_one({"task_id": task_id, "level": "G"})
    except Exception as e:
        print  e

        # bind_obj.set_error("cannot find annotation_query main table", code="51008159")

    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['exp_detail']
    level = bind_obj.sheet.option('level')
    if level[0].upper() == 'T':
        query_dict = dict(exp_id=annot_main_id,)
        result_dict = dict( _id=0, gene_name=1, gene_id=1, transcript_id=1)
        result = annot_detail.find(query_dict, result_dict)
        gene2name = pd.DataFrame(list(result))
        gene2name.set_index('transcript_id', inplace=True)
    else:
        query_dict = dict(exp_id=annot_main_id)
        result_dict = dict( _id=0, gene_name=1, gene_id=1)
        result = annot_detail.find(query_dict, result_dict)
        gene2name = pd.DataFrame(list(result))
        gene2name.set_index('gene_id', inplace=True)
    gene2name = gene2name.loc[list(target_seqs), :]
    gene2name.reset_index(inplace=True)
    output = os.path.join(dir_path, "seq_id2gene_name.txt")
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    gene2name = pd.read_table(output, header=0)
    gene2name = gene2name.fillna(method="pad", axis=1)
    gene2name.to_csv(output, sep='\t', header=True, index=False)
    return exp_matrix+';'+output


# added  by gdq for wgcna
def export_wgcna_relate_input(data, option_name, dir_path, bind_obj=None):
    module_id = data.strip()
    # export eigengenes
    eigengenes = db['wgcna_module_eigengenes_detail']
    eigengenes_found = eigengenes.find({"module_id": ObjectId(module_id)}, {"_id": 0, "module_id":0})
    eigengenes_pd = pd.DataFrame(list(eigengenes_found))
    eigengenes_pd.set_index("module", inplace=True)
    eigengenes_path = os.path.join(dir_path, "module_eigengenes.xls")
    eigengenes_pd.to_csv(eigengenes_path, sep='\t', header=True, index=True)
    # export exp matrix
    prepare_id = db['wgcna_module'].find_one({"main_id": ObjectId(module_id)})['wgcna_prepare_id']
    exp_matrix = db['wgcna_prepare'].find_one({"main_id": ObjectId(prepare_id)})['exp_matrix']
    # export each gene's module and gene_id/gene_name info
    level = bind_obj.sheet.option('level')
    membership = db['wgcna_module_membership_detail']
    query_dict = dict(module_id=ObjectId(module_id))
    return_dict = dict(_id=0, seq_id=1, gene_name=1, module=1, kme=1, block_id=1)
    if level == 'T':
        return_dict.update({"gene_id": 1})
    result = membership.find(query_dict, return_dict)
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("seq_id", inplace=True)
    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    return exp_matrix + ';' + eigengenes_path + ";" + gene_annot


def get_after_qc_bam_path(data, option_name, dir_path, bind_obj=None):
    db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
    conn = db['specimen']
    result=conn.find({"task_id": data, "about_qc": "after"})
    sorted_result = sorted(result, key=lambda k: k['_id'])
    output = os.path.join(dir_path, option_name)
    with open(output,'w') as f:
        for i in sorted_result:
            f.write(i['bam_path'] + "\n")
    return output

#---------------------从蛋白复制过来的------------------------

def export_compare_exp_fc(data, option_name, dir_path, bind_obj=None):
    '''
    导出差异分析accession_id 和 fc
    '''

    diff_id = bind_obj.sheet.option('diff_id')
    compare_group = bind_obj.sheet.option('compare_group')
    target_cols = OrderedDict(seq_id=1, log2fc=1, _id=0)

    bind_obj.logger.debug("导出表达参数 {}".format(target_cols))
    conn = db['diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id),"compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    output = os.path.join(dir_path, option_name)
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export expression matrix')
    return output

def export_go_enrich_matrix(data, option_name, dir_path, bind_obj=None):

    go_enrich_id = bind_obj.sheet.option('enrich_id')
    go_type = bind_obj.sheet.option('go_type')

    target_cols = OrderedDict(go_id=1, go_type=1, discription=1, p_corrected=1, p_uncorrected=1, seq_list=1, depth=1, _id=0)
    bind_obj.logger.debug("导出GO富集表")
    # get geneset
    conn = db['geneset_go_enrich_detail']
    if go_type == 'ALL' or go_type == 'All':
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id)}, target_cols)
    else:
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id), "go_type":go_type}, target_cols)
    go_enrich_matrix = pd.DataFrame(list(go_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in go_enrich_matrix['seq_list']]
    go_enrich_matrix['seq_list'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    go_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export go_enrich matrix')
    return output

def export_kegg_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    导出kegg富集表格
    '''
    kegg_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, id=1, term=1, pvalue=1, corrected_pvalue=1, seq_list=1, kegg_type=1)
    bind_obj.logger.debug("导出KEGG富集表")
    # get geneset
    conn = db['geneset_kegg_enrich_detail']
    bind_obj.logger.debug("导出以下列{}".format(target_cols))
    bind_obj.logger.debug("导出以下列{}".format(kegg_enrich_id))
    kegg_enrich_records = conn.find({"kegg_enrich_id": ObjectId(kegg_enrich_id)}, target_cols)
    kegg_enrich_matrix = pd.DataFrame(list(kegg_enrich_records))
    seq_list2 = [i.replace(';', '|') for i in kegg_enrich_matrix['seq_list']]
    kegg_enrich_matrix['seq_list'] = seq_list2
    # exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    kegg_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export kegg_enrich matrix')
    return output

def export_gene_list_ppi(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "%s.txt" % option_name)
    bind_obj.logger.debug("正在导出蛋白集")
    collection = db['geneset_detail']
    main_collection = db['geneset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:%s在geneset中未找到！", variables=(ObjectId(data)), code="51008160")
    results = collection.find_one({"geneset_id": ObjectId(data)})["seq_list"]
    with open(gene_list_path, "wb") as f:
        f.write("accession_id" + "\n")
        for result in results:
            f.write(result + "\n")
    bind_obj.logger.debug("蛋白集导出成功！")
    return gene_list_path

def export_geneset_genes_name(data, option_name, dir_path, bind_obj=None):
    exp_id, geneset_id = data.split(";")
    task_id = bind_obj.sheet.option('task_id')
    bind_obj.logger.info("sheet is {}".format(bind_obj.sheet.option('level')))
    try:
        level = bind_obj.sheet.option('level')

        if level == "G":
            id_name = "gene_id"
        else:
            id_name = "transcript_id"
        # bind_obj.logger.info("level is {}".format(level))
    except:
        level = "G"
        id_name = "gene_id"
    if "all" not in geneset_id.lower():
        conn = db['geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        if not geneset_records:
            bind_obj.set_error('geneset not found by query: %s', variables=(geneset_id), code="51008154")
        seq_id_set = set(geneset_records['seq_list'])
    else:
        seq_id_set = set()

    col = db["exp"]
    main_id = col.find_one({"task_id": task_id, "level": level})['main_id']

    bind_obj.logger.info("level is {}".format(level))

    col = db["exp_detail"]
    select_dict = {'gene_id': 1, 'gene_name': 1, '_id': 0}
    if level == "T":
        select_dict = {'transcript_id': 1, 'gene_name': 1, '_id': 0}
        if "wgcna" in bind_obj.sheet.name:
            select_dict = {'transcript_id': 1, 'gene_id': 1, 'gene_name': 1, '_id': 0}
    bind_obj.logger.info("select_dict is {}".format(select_dict))
    results = col.find({"exp_id": ObjectId(main_id)}, select_dict)
    out_file = os.path.join(dir_path, 'geneset.name')
    if level == "T" and "wgcna" in bind_obj.sheet.name:
        line_demo = '{transcript_id}\t{gene_id}\t{gene_name}\n'
    else:
        line_demo = '{gene_id}\t{gene_name}\n'
    with open(out_file, 'w') as out_handler:
        if level == "T" and "wgcna" in bind_obj.sheet.name:
            out_handler.write('transcript_id\tgene_id\tgene_name\n')
        else:
            out_handler.write('gene_id\tgene_name\n')
        for dic in results:
            if level == "T" and "wgcna" in bind_obj.sheet.name:
                out_handler.write(line_demo.format(transcript_id=dic["transcript_id"], gene_id=dic.get("gene_id", ""), gene_name=dic.get("gene_name", dic[id_name])))

            else:
                if len(seq_id_set) == 0:
                    out_handler.write(line_demo.format(gene_id=dic[id_name], gene_name=dic.get("gene_name", dic[id_name])))
                else:
                    if dic[id_name] in seq_id_set:
                        out_handler.write(line_demo.format(gene_id=dic[id_name], gene_name=dic.get("gene_name", dic[id_name])))
    return out_file


def export_genes_detail(data, option_name, dir_path, bind_obj=None):
    task_id, type_ = data['task_id'], data['type']
    level = bind_obj.sheet.option("level")
    col = db["exp"]
    main_id = col.find_one({"task_id": task_id, "level": level})['main_id']

    col = db["exp_detail"]
    if level == "G" :
        choose_fields = {'gene_id': 1, 'gene_name': 1, 'description': 1, '_id': 0}
    else:
        choose_fields = {'transcript_id': 1, 'gene_name': 1, 'description': 1, '_id': 0}

    results = col.find({"exp_id": ObjectId(main_id)}, choose_fields)
    # Probe Set ID    Gene Symbol     Gene Title
    out_list = []
    line_demo = '{gene_id}\t{gene_name}\t{description}\n'
    seq_id_set = set()
    out_file = os.path.join(dir_path, 'gsea.chip')
    try:
        geneset_source = bind_obj.sheet.option('geneset_source')
    except:
        geneset_source = ""


    with open(out_file, 'w') as out_handler:
        out_handler.write('gene_id\tgene_name\tdescription\n')
        for dic in results:
            if "gene_name" in dic:
                pass
            else:
                continue

            if geneset_source == 'msigdb' and dic["gene_name"] in [None, "None", "", "-"]:
                continue
            if "description" in dic:
                if dic["description"] in [None, "None", ""]:
                    dic["description"] = "-"
            else:
                dic["description"] = "-"
            if level == "T":
                gene_id = dic['transcript_id']
            else:
                gene_id = dic['gene_id']
            if gene_id in seq_id_set:
                continue
            seq_id_set.add(gene_id)
            # bind_obj.logger.info("data is {}".format(dic))
            if level == "G":
                out_handler.write(line_demo.format(gene_id=dic["gene_id"], gene_name=dic["gene_name"], description=dic["description"]))
            if level == "T":
                out_handler.write(line_demo.format(gene_id=dic["transcript_id"], gene_name=dic["gene_name"], description=dic["description"]))

    return out_file

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

# added for transcription factor analysis  --------------------------------------------
def get_all_pep_seq(data, option_name, dir_path, bind_obj=None):
    pep_db_path = data.strip()
    result_path = os.path.join(dir_path, "all_pep.fa")
    if os.path.exists(pep_db_path):
        pass
    elif re.match(r'^\w+://\S+/.+$', pep_db_path) or re.match(r'/mnt/ilustre', pep_db_path):
        transfer = MultiFileTransfer()
        transfer.add_download(pep_db_path, bind_obj.work_dir + "/")
        transfer.perform()
        #download_s3_file(pep_db_path, os.path.join(bind_obj.work_dir, "pep.db"))
        base_name = os.path.basename(pep_db_path)
        pep_db_path = os.path.join(bind_obj.work_dir, base_name)
    bind_obj.logger.debug("文件为 {}".format(pep_db_path))
    import sqlite3
    conn = sqlite3.connect(pep_db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT transcript_id,pep_seq FROM sub_info")
    with open(result_path, 'w') as fw:
        # pep_id 就是转录本的id
        for pep_id, pep_seq in cursor.fetchall():
            fw.write('>{}\n{}\n'.format(pep_id, pep_seq))
    return result_path


def get_gene_detail(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip()
    exp_table = db['exp']
    exp_main = exp_table.find_one({"task_id": task_id, "level": "T", 'status': 'end'})
    if not exp_main:
        bind_obj.set_error("Not Found in exp by query %s", variables=(task_id))
    if "main_id" not in exp_main:
        exp_main_id = exp_main['_id']
    else:
        exp_main_id = exp_main['main_id']
    exp_detail = db['exp_detail']
    query_dict = dict(exp_id=exp_main_id, )
    result_dict = dict(_id=0, gene_name=1, gene_id=1, description=1, transcript_id=1)
    result = exp_detail.find(query_dict, result_dict)
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("transcript_id", inplace=True)
    #result_pd = result_pd.loc[:, ["gene_id", "gene_name", "description"]]
    #result_pd.columns = ["gene_id", "gene_name", "gene_desc"]
    result_pd = result_pd.loc[:, ["gene_id", "gene_name"]]
    result_pd.columns = ["gene_id", "gene_name"]
    result_pd["gene_name"][result_pd['gene_name'] == "-"] = None
    result_pd = result_pd.fillna(method="pad", axis=1)
    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    os.system(r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot

def get_gene_detail_new(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip().split("|")[0]
    geneset_type = data.strip().split("|")[1]
    annot_table = db['annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", 'status': 'end'})
    if not annot_main:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    if not annot_main:
        bind_obj.set_error("Not Found in annotation_query by query %s", variables=(task_id), code="51008162")
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['annotation_query_detail']
    query_dict = dict(query_id=annot_main_id, )
    if geneset_type.lower() == "t":
        result_dict = dict(_id=0, gene_name=1, gene_id=1, description=1, transcript_id=1)
        result = annot_detail.find(query_dict, result_dict)
        result_pd = pd.DataFrame(list(result))
        result_pd.set_index("transcript_id", inplace=True)
        result_pd = result_pd.loc[:, ["gene_id", "gene_name", "description"]]
        result_pd.columns = ["gene_id", "gene_name", "gene_desc"]
        result_pd["gene_name"][result_pd['gene_name'] == "-"] = None
        result_pd = result_pd.fillna(method="pad", axis=1)
    else:
        result_dict = dict(_id=0, gene_name=1, gene_id=1, description=1)
        result = annot_detail.find(query_dict, result_dict)
        result_pd = pd.DataFrame(list(result))
        result_pd.set_index("gene_id", inplace=True)
        result_pd = result_pd.loc[:, ["gene_name", "description"]]
        result_pd.columns = ["gene_name", "gene_desc"]
        result_pd["gene_name"][result_pd['gene_name'] == "-"] = None
        result_pd = result_pd.fillna(method="pad", axis=1)
    gene_annot = os.path.join(dir_path, "seq_annot.xls")
    result_pd.to_csv(gene_annot, sep='\t', header=True, index=True)
    os.system(r"sed -i 's/%2B/+/g;s/%2F/\//g;s/%2C/,/g;s/%3A/:/g;s/%3B/;/g;s/%3D/=/g;s/%3F/?/g;s/%20/ /g;s/%25/%/g;s/%3C/</g;s/%3E/>/g;s/%5B/[/g;s/%5D/]/g;s/%7B/{/g;s/%7D/}/g' " + gene_annot)
    return gene_annot

# added by gdq for tfbs predict
def export_geneid2tfid_file(data, option_name, dir_path, bind_obj=None):
    conn = db['tf_predict_detail']
    result = conn.find({"tf_predict_id": ObjectId(data.strip())}, {"_id": 0, 'gene_id': 1, "blast_hit": 1})
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("gene_id", inplace=True)
    result_file = os.path.join(dir_path, "geneid2tfid.txt")
    result_pd.to_csv(result_file, sep='\t', header=True, index=True)
    return result_file


def export_predict_result(data, option_name, dir_path, bind_obj=None):
    conn = db['tfbs_predict_detail']
    target_cols = {
        "_id":0, "tf_geneid":1, "gene_name_tf":1, "target_id":1, "gene_name_target":1,
        "corr":1, "corr_pvalue":1, "corr_padjust":1, "p-value":1, "q-value":1
    }
    result = conn.find({"tfbs_predict_id": ObjectId(data.strip()),}, target_cols,)
    result_pd = pd.DataFrame(list(result))
    result_file = os.path.join(dir_path, "tfbs_predict.xls")
    result_pd.to_csv(result_file, sep='\t', header=True, index=False)
    return result_file

# ---------------------------------------------------------------------------------------------

def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} = {}'.format(k, v))

def export_stem_matrix(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    collection = db['exp_detail']
    cursor = collection.find({'exp_id': ObjectId(data)})
    output = os.path.join(dir_path, 'matrix.tsv')
    df = pd.DataFrame(list(cursor))
    df.to_csv(output, sep='\t', index=False)
    return output

def export_stem_geneset(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    output = os.path.join(dir_path, 'geneset.list')
    if data in ['All', 'RefAll']:
        collection = db['exp_detail']
        cursor = collection.find({'exp_id': ObjectId(bind_obj.sheet.option('matrix'))})
        if data == 'All':
            dct = {'seq_list': list({document['seq_id'] for document in cursor})}
        elif data == 'RefAll':
            dct = {'seq_list': list({document['seq_id'] for document in cursor if not document['is_new']})}
    else:
        collection = db['geneset_detail']
        dct = collection.find_one({'geneset_id': ObjectId(data)})
    open(output, 'w').writelines('{}\n'.format(seq) for seq in dct['seq_list'])
    return output

def export_stem_group(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    collection = db['specimen_group']
    dct = collection.find_one({'main_id': ObjectId(data)})
    group_list = dct['category_names']
    samples_list = dct['specimen_names']
    lines = ['#sample\tgroup\n']
    for samples, group in zip(samples_list, group_list):
        lines.extend('{}\t{}\n'.format(sample, group) for sample in samples)
    else:
        output = os.path.join(dir_path, 'group.txt')
        open(output, 'w').writelines(lines)
        return output

def export_masigpro_geneset(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    output = os.path.join(dir_path, 'geneset.list')
    if data in ['All', 'RefAll']:
        collection = db['exp_detail']
        cursor = collection.find({'exp_id': ObjectId(bind_obj.sheet.option('matrix'))})
        if data == 'All':
            dct = {'seq_list': list({document['seq_id'] for document in cursor})}
        elif data == 'RefAll':
            dct = {'seq_list': list({document['seq_id'] for document in cursor if not document['is_new']})}
    else:
        collection = db['geneset_detail']
        dct = collection.find_one({'geneset_id': ObjectId(data)})
    open(output, 'w').writelines('{}\n'.format(seq) for seq in dct['seq_list'])
    return output

def export_masigpro_matrix(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    collection = db['exp_detail']
    cursor = collection.find({'exp_id': ObjectId(data)})
    output = os.path.join(dir_path, 'matrix.tsv')
    df = pd.DataFrame(list(cursor))
    df.to_csv(output, sep='\t', index=False)
    return output


def export_exp_matrix(data, option_name, dir_path, bind_obj=None):
    exp_id = data
    group_dict = json.loads(bind_obj.sheet.option('group'), object_pairs_hook=OrderedDict)
    geneset_id = bind_obj.sheet.option('geneset_id')
    if geneset_id.lower() in ["all", "refall"]:
        pass
    else:
        seq_list = db['geneset_detail'].find_one({'geneset_id': ObjectId(geneset_id)})['seq_list']

    samples = list()
    for i in group_dict.values():
        samples.extend(i)

    col = db["exp_detail"]
    if geneset_id.lower() == "refall":
        results = col.find({"exp_id": ObjectId(exp_id), 'is_new': False}, {'_id': 0, 'exp_id': 0, 'is_new': 0})
    else:
        results = col.find({"exp_id": ObjectId(exp_id)}, {'_id': 0, 'exp_id': 0, 'is_new': 0})


    df = pd.DataFrame(list(results))
    df.set_index('gene_id', inplace=True)
    if geneset_id.lower() in ["all", "refall"]:
        df = df[samples]
    else:
        df = df[samples].loc[df.index & set(seq_list)]
    df.index.name = 'seq_id'
    out_file = os.path.join(dir_path, 'gsea.txt')
    df.to_csv(out_file, sep='\t')

    return out_file


def export_anno_genesets(data, option_name, dir_path, bind_obj=None):
    # {'c1': data.c1, 'c2': data.c2, 'c3': data.c3, 'task_id': data.task_id,
    #                                   'type': 'origin'}
    bind_obj.logger.info("data is {}".format(data))
    task_id = bind_obj.sheet.option('task_id')
    anno_level = bind_obj.sheet.option('level')
    if (data.startswith("GO")):
        anno_type, subtype, genesets = data.split(";")
    else:
        anno_type, subtype, subtype2, genesets = data.split(";")
    go = {
        'MF': 'molecular_function',
        'CC': 'cellular_component',
        'BP': 'biological_process'
    }
    kegg = {
        "Environmental Information Processing",
        "Cellular Processes",
        "Organismal Systems",
        "Human Diseases",
        "Metabolism",
        "Genetic Information Processing"
    }
    if anno_type.lower() == 'go':
        main_id = db['annotation_go'].find_one({'task_id': task_id})['main_id']
        query_dic = {'go_id': main_id, 'anno_type': anno_level}
        '''
        if subtype.upper() in go:
            query_dic['goterm'] = go[subtype.upper()]
        '''
        results = db['annotation_go_list'].find(query_dic, {'gene_id': 1, 'gos_list': 1, '_id': 0})
        '''
        go_list = None
        if genesets.strip():
            go_list = {i.strip() for i in genesets.strip().split(',')}
        '''
        out_file = os.path.join(dir_path, 'go.list')
        with open(out_file, 'w') as out_handler:
            out_handler.write(
                ''.join(['{}\t{}\n'.format(dic['gene_id'], dic['gos_list']) for dic in results])
            )
        return out_file
    else:
        main_id = db['annotation_kegg'].find_one({'task_id': task_id})['main_id']
        query_dic = {'kegg_id': main_id, 'anno_type': anno_level}
        if subtype == "defined":
            pass
        else:
            query_dic['first_category'] = subtype
            if subtype2 not in ['all', 'All', ""]:
                if len(subtype2.split(",")) > 1:
                    subtype2_par = subtype2.split(",")
                    subtype2_list = [s.replace("_", ",") for s in subtype2_par]
                    query_dic['second_category'] = {"$in": subtype2_list}

                else:
                    query_dic['second_category'] = subtype2.replace("_", ",")


        query_dic['seq_type'] = "all"

        bind_obj.logger.info("find {} {}".format('annotation_kegg_level', query_dic))
        results = db['annotation_kegg_level'].find(query_dic, {'seq_list': 1,'pathway_definition': 1, 'pathway_id': 1, '_id': 0})
        kegg_list = None
        if genesets not in ["all", ""]:
            kegg_list = {i.strip() for i in genesets.strip().split(',')}
        out_file = os.path.join(dir_path, 'kegg.gmt')
        with open(out_file, 'w') as out_handler:
            if kegg_list:
                out_handler.write(
                    ''.join('{}\t{}\t{}\n'.format(
                        dic['pathway_id'], dic['pathway_definition'], dic['seq_list'].replace(';', '\t')) for dic in results if dic['pathway_id'] in kegg_list)
                )
            else:
                out_handler.write(
                    ''.join(['{}\t{}\t{}\n'.format(
                        dic['pathway_id'], dic['pathway_definition'], dic['seq_list'].replace(';', '\t')) for dic in results])
                )
        return out_file


def export_msigdb_genesets(data, option_name, dir_path, bind_obj=None):
    # {'c1': data.c1, 'c2': data.c2, 'c3': data.c3}
    bind_obj.logger.info("data is {}".format(data))
    extr_dic = {}
    species= data.split(";")[0]
    extr_dic['organism'] = species
    c1 = data.split(";")[1]
    extr_dic['c1'] = c1
    c2 = data.split(";")[2]
    if c2.count('|') == 1:
        c2, tmp_c3 = c2.split('|')
        extr_dic['c2'] = c2
        extr_dic['c3'] = tmp_c3
    else:
        extr_dic['c2'] = c2

    results = db['msigdb'].find(extr_dic, {'_id': 0})
    ref_link = re.compile('<[^>]*>')
    res_dict = {dic['name']: '{}\t{}\t{}\n'.format(
        dic["name"], re.sub(ref_link, '', dic['brief_description']), dic["gene_symbols"].replace(" ", "\t")) for dic in results if 'name' in dic}
    c3 = data.split(";")[3]
    out_file = os.path.join(dir_path, 'gsea.gmt')
    with open(out_file, 'w') as out_handler:
        if c3 == 'all' or c3 =='':
            out_handler.write(''.join(res_dict.values()))
        else:
            keys = {i.strip() for i in c3.split(',')}
            out_handler.write(''.join(v for k, v in res_dict.items() if k in keys))
    return out_file

# ---------------------------------------------------------------------------------------------
