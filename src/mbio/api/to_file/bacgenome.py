# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
import os
import re
import json
from types import StringTypes
from biocluster.config import Config
from bson.objectid import ObjectId
from collections import defaultdict
from mainapp.models.mongo.bacgenome import Bacgenome
from mainapp.controllers.project.bacgenome_controller import BacgenomeController

client = Config().get_mongo_client(mtype="bacgenome")
db = client[Config().get_mongo_dbname("bacgenome")]


def export_gene_by_geneid(data, option_name, dir_path, bind_obj=None):
    """
    根据gene_id/样品名生成基因序列文件
    :return:
    """
    file_path = os.path.join(dir_path, "%s.fasta" % option_name)
    collection_main = db['gene_predict']
    main = collection_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
    main_id = main['_id']
    collection_seq = db['gene_predict_seq']
    doc = collection_seq.find_one(
        {'predict_id': main_id, 'specimen_id': bind_obj.sheet.option("specimen_id"), 'gene_id': data})
    with open(file_path, 'w') as f:
        f.write(doc['seq_info'] + '\n' + doc['seq_fnn'] + '\n')
    return file_path

def export_gene_faa_by_geneid(data, option_name, dir_path, bind_obj=None):     #guanqing.zou 20180912
    """
    根据gene_id/样品名生成基因faa序列文件
    :return:
    """
    file_path = os.path.join(dir_path, "%s.fasta" % option_name)
    collection_main = db['gene_predict']
    main = collection_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
    main_id = main['_id']
    collection_seq = db['gene_predict_seq']
    doc = collection_seq.find_one(
        {'predict_id': main_id, 'specimen_id': bind_obj.sheet.option("specimen_id"), 'gene_id': data})
    with open(file_path, 'w') as f:
        f.write(doc['seq_info'] + '\n' + doc['seq_faa'] + '\n')
    return file_path

def export_prephage_circos(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "%s.txt" % option_name)
    assemble = db['assemble']
    main_assemble = assemble.find_one({'task_id': bind_obj.sheet.option("task_id")})
    if main_assemble['data_type'] == 'annotation':
        is_anno_pipline = True
    else:
        is_anno_pipline = False
    analysis_type = main_assemble['analysis_type']

    assemble_id = main_assemble['_id']
    seq_collection = db['assemble_seq']

    def get_pos(num, pos):
        sum = 0
        for i in range(1, num):
            if is_anno_pipline:
                one_seq = seq_collection.find_one(
                    {'assemble_id': assemble_id, 'specimen_id': bind_obj.sheet.option("specimen_id"),
                    'seq_sort_id': 'seq' + str(i)})
            else:
                one_seq = seq_collection.find_one(
                    {'assemble_id': assemble_id, 'specimen_id': bind_obj.sheet.option("specimen_id"),
                    'seq_id': 'Scaffold' + str(i)})
            sum += one_seq['len']
        new_pos = sum + pos
        return new_pos

    type = data.split(",")
    with open(file_path, 'w') as f:
        if "prephage" in type:
            collection_main = db['prephage']
            main = collection_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
            main_id = main['_id']
            collection_detail = db['prephage_detail']

            #if bind_obj.sheet.option("location").startswith('p') or bind_obj.sheet.option("location").startswith('Chr'):
            if analysis_type == 'complete':
                doc = collection_detail.find(
                    {'prephage_id': main_id, 'specimen_id': bind_obj.sheet.option("specimen_id"),
                     'location': bind_obj.sheet.option("location")})
                for one in doc:
                    start = int(one['ph_start'])
                    end = int(one['ph_end'])
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join([one['location'], str(int(s)), str(int(e)), 'fill_color=(123,104,238)']) + '\n')
            else:
                doc = collection_detail.find(
                    {'prephage_id': main_id, 'specimen_id': bind_obj.sheet.option("specimen_id")})
                for one in doc:
                    if is_anno_pipline:
                        scaffold_num = int(one['location'][3:])  #seq1
                    else:
                        scaffold_num = int(one['location'][8:])  #scaffold1
                    start = get_pos(scaffold_num, int(one['ph_start']))
                    end = get_pos(scaffold_num, int(one['ph_end']))
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join(['Scaffold', str(int(s)), str(int(e)), 'fill_color=(123,104,238)']) + '\n')
        if "gi" in type:
            collection_main = db['island']
            main = collection_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
            main_id = main['_id']
            collection_detail = db['island_detail']
            #if bind_obj.sheet.option("location").startswith('p') or bind_obj.sheet.option("location").startswith('Chr'):
            if analysis_type == 'complete':
                doc = collection_detail.find(
                    {'island_id': main_id, 'specimen_id': bind_obj.sheet.option("specimen_id"),
                     'location': bind_obj.sheet.option("location")})
                for one in doc:
                    start = int(one['island_start'])
                    end = int(one['island_end'])
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join([one['location'], str(int(s)), str(int(e)), 'fill_color=(205,205,0)']) + '\n')
            else:
                doc = collection_detail.find(
                    {'island_id': main_id, 'specimen_id': bind_obj.sheet.option("specimen_id")})
                for one in doc:
                    if is_anno_pipline:
                        scaffold_num = int(one['location'][3:])  #seq1
                    else:
                        scaffold_num = int(one['location'][8:])
                    start = get_pos(scaffold_num, int(one['island_start']))
                    end = get_pos(scaffold_num, int(one['island_end']))
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join(['Scaffold', str(int(s)), str(int(e)), 'fill_color=(205,205,0)']) + '\n')
        if "ncrna" in type:
            rrna_predict_main = db['rrna_predict']
            rrna_predict_main = rrna_predict_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
            rrna_id = rrna_predict_main['_id']
            rrna_detail = db['rrna_predict_detail']
            trna_predict_main = db['trna_predict']
            trna_predict_main = trna_predict_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
            trna_id = trna_predict_main['_id']
            trna_detail = db['trna_predict_detail']
            #if bind_obj.sheet.option("location").startswith('p') or bind_obj.sheet.option("location").startswith('Chr'):
            if analysis_type == 'complete':
                doc = trna_detail.find(
                    {'predict_id': trna_id, 'specimen_id': bind_obj.sheet.option("specimen_id"),
                     'location': bind_obj.sheet.option("location")})
                for one in doc:
                    start = int(one['start'])
                    end = int(one['end'])
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join([one['location'], str(int(s)), str(int(e)), 'fill_color=(255,0,0)']) + '\n')
                docr = rrna_detail.find(
                    {'predict_id': rrna_id, 'specimen_id': bind_obj.sheet.option("specimen_id"),
                     'location': bind_obj.sheet.option("location")})
                for one in docr:
                    start = int(one['start'])
                    end = int(one['end'])
                    s = min(start, end)
                    e = max(start, end)
                    if one['type'] == '5S_rRNA':
                        color = '(0,128,0)'
                    elif one['type'] == '16S_rRNA':
                        color = '(0,0,255)'
                    elif one['type'] == '23S_rRNA':
                        color = '(0,191,255)'
                    else:
                        raise Exception('wrong rRNA type')
                    f.write(' '.join([one['location'], str(int(s)), str(int(e)), 'fill_color=' + color]) + '\n')

            else:
                doc = trna_detail.find(
                    {'predict_id': trna_id, 'specimen_id': bind_obj.sheet.option("specimen_id")})
                for one in doc:
                    if is_anno_pipline:
                        scaffold_num = int(one['location'][3:])  #seq1
                    else:
                        scaffold_num = int(one['location'][8:])
                    start = get_pos(scaffold_num, int(one['start']))
                    end = get_pos(scaffold_num, int(one['end']))
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join(['Scaffold', str(int(s)), str(int(e)), 'fill_color=(255,0,0)']) + '\n')
                docr = rrna_detail.find(
                    {'predict_id': rrna_id, 'specimen_id': bind_obj.sheet.option("specimen_id")})
                for one in docr:
                    if is_anno_pipline:
                        scaffold_num = int(one['location'][3:])  #seq1
                    else:
                        scaffold_num = int(one['location'][8:])
                    start = get_pos(scaffold_num, int(one['start']))
                    end = get_pos(scaffold_num, int(one['end']))
                    s = min(start, end)
                    e = max(start, end)
                    if one['type'] == '5S_rRNA':
                        color = '(0,128,0)'
                    elif one['type'] == '16S_rRNA':
                        color = '(0,0,255)'
                    elif one['type'] == '23S_rRNA':
                        color = '(0,191,255)'
                    else:
                        raise Exception('wrong rRNA type')
                    f.write(
                        'Scaffold' + ' ' + str(int(s)) + ' ' + str(int(e)) + ' ' + 'fill_color=' + color + '\n')
        if "is" in type:
            collection_main = db['is']
            main = collection_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
            main_id = main['_id']
            collection_detail = db['is_detail']

            # if bind_obj.sheet.option("location").startswith('p') or bind_obj.sheet.option("location").startswith('Chr'):
            if analysis_type == 'complete':
                doc = collection_detail.find(
                    {'is_id': main_id, 'sample': bind_obj.sheet.option("specimen_id"),
                     'location': bind_obj.sheet.option("location")})
                for one in doc:
                    start = int(one['start'])
                    end = int(one['end'])
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join([one['location'], str(int(s)), str(int(e)), 'fill_color=(28,28,28)']) + '\n')
            else:
                doc = collection_detail.find(
                    {'is_id': main_id, 'sample': bind_obj.sheet.option("specimen_id")})
                for one in doc:
                    if is_anno_pipline:
                        scaffold_num = int(one['location'][3:])  # seq1
                    else:
                        scaffold_num = int(one['location'][8:])  # scaffold1
                    start = get_pos(scaffold_num, int(one['start']))
                    end = get_pos(scaffold_num, int(one['end']))
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join(['Scaffold', str(int(s)), str(int(e)), 'fill_color=(28,28,28)']) + '\n')
        if "integron" in type:
            collection_main = db['integron']
            main = collection_main.find_one({'task_id': bind_obj.sheet.option("task_id")})
            main_id = main['_id']
            collection_detail = db['integron_detail']

            # if bind_obj.sheet.option("location").startswith('p') or bind_obj.sheet.option("location").startswith('Chr'):
            if analysis_type == 'complete':
                doc = collection_detail.find(
                    {'integron_id': main_id, 'sample': bind_obj.sheet.option("specimen_id"),
                     'location': bind_obj.sheet.option("location")})
                for one in doc:
                    start = int(one['start'])
                    end = int(one['end'])
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join([one['location'], str(int(s)), str(int(e)), 'fill_color=(255,127,0)']) + '\n')
            else:
                doc = collection_detail.find(
                    {'integron_id': main_id, 'sample': bind_obj.sheet.option("specimen_id")})
                for one in doc:
                    if is_anno_pipline:
                        scaffold_num = int(one['location'][3:])  # seq1
                    else:
                        scaffold_num = int(one['location'][8:])  # scaffold1
                    start = get_pos(scaffold_num, int(one['start']))
                    end = get_pos(scaffold_num, int(one['end']))
                    s = min(start, end)
                    e = max(start, end)
                    f.write(' '.join(['Scaffold', str(int(s)), str(int(e)), 'fill_color=(255,127,0)']) + '\n')
        os.system('sort ' + file_path)
    return file_path

def export_def_circle(data, option_name, dir_path, bind_obj=None):   #zouguanqing 20190410
    file_path = os.path.join(dir_path, "%s.txt" % option_name)
    assemble = db['assemble']
    main_assemble = assemble.find_one({'task_id': bind_obj.sheet.option("task_id")})
    assemble_id = main_assemble['_id']
    seq_collection = db['assemble_seq']

    if main_assemble['data_type'] == 'annotation':
        is_anno_pipline = True
    else:
        is_anno_pipline = False
    analysis_type = main_assemble['analysis_type']

    def get_pos(num, pos):
        sum = 0
        for i in range(1, num):
            if is_anno_pipline:
                one_seq = seq_collection.find_one(
                    {'assemble_id': assemble_id, 'specimen_id': bind_obj.sheet.option("specimen_id"),
                    'seq_sort_id': 'seq' + str(i)})
            else:
                one_seq = seq_collection.find_one(
                    {'assemble_id': assemble_id, 'specimen_id': bind_obj.sheet.option("specimen_id"),
                    'seq_id': 'Scaffold' + str(i)})
            if not one_seq:
                bind_obj.logger.error('assemble_seq database not find info')
            sum += one_seq['len']
        new_pos = sum + pos
        return new_pos
    color ='(255,0,255)'
    temp = list()
    # if 'Scaffold' in data:
    #     seq_name = 'Scaffold'
    # elif 'Chromosome' in data:
    #     seq_name = 'Chr'
    # else:
    seq_name = data.split(';')[0].split(',')[0]
    with open(file_path, 'w') as f:
        all_data = data.split(';')
        for ds in all_data:
            d = ds.split(',')
            if analysis_type != 'complete':
                if is_anno_pipline:
                    num = int(d[0][3:])
                else:
                    num = int(d[0][8:])
                new_pos1 = int(get_pos(num,int(d[1])))
                new_pos2 = int(get_pos(num,int(d[2])))

            else:
                new_pos1 = int(d[1])
                new_pos2 = int(d[2])

            if new_pos1 > new_pos2:
                new_pos1, new_pos2 = new_pos2, new_pos1
            temp.append([new_pos1,new_pos2])
        stemp = sorted(temp, key=lambda a:a[0])
        for i in stemp:
            f.write(seq_name + ' ' + str(i[0]) + ' ' + str(i[1]) + ' ' + 'fill_color=' + color + '\n')
    return file_path

# 创建基因集文件
def export_geneset(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    geneset_id = ObjectId(data)
    gene_predict_tool = db['gene_predict_tool']
    main_gene_predict_tool = gene_predict_tool.find_one({'_id': geneset_id})
    geneset = main_gene_predict_tool["geneset"].split(",")
    specimen_id = main_gene_predict_tool["specimen_id"]

    file_path = os.path.join(dir_path, "%s.txt" % option_name)
    gene_predict = db['gene_predict']
    main_gene_predict = gene_predict.find_one({'task_id': bind_obj.sheet.option("task_id")})
    gene_predict_id = main_gene_predict['_id']
    gene_predict_seq = db['gene_predict_seq']
    all_gene = gene_predict_seq.find({'predict_id': gene_predict_id,'specimen_id':specimen_id})
    with open(file_path,"w") as t:
        t.write("gene_id\tgene_name\tseq\n")
        for gene in all_gene:
            if gene["gene_id"] in geneset:
                t.write(gene["gene_id"]+"\t"+gene["gene_id"]+"\t"+gene["seq_fnn"]+"\n")
    return file_path
