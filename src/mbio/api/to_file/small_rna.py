# -*- coding: utf-8 -*-
# __author__ = "sanger"

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
import sys
import numpy as np

project_type = 'small_rna'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

def export_kegg_pdf(data, option_name, dir_path, bind_obj=None):
    fs = gridfs.GridFS(db)
    annotation_collection = db["sg_annotation_kegg"]
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    main_id = annotation_collection.find_one({"task_id":task_id, "type": bind_obj.sheet.option("type")})["main_id"]
    kegg_level_collection = db["sg_annotation_kegg_level"]
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

# ---------------------------------------- geneset related qdh -----------------------------------------

def get_geneset_detail(data, bind_obj):
    geneset_collection = db['sg_geneset']
    genesets = dict()
    names = list()
    task_id = str()
    geneset_type = bind_obj.sheet.option('geneset_type')
    for geneset_id in data.split(','):
        geneset_result = geneset_collection.find_one({'main_id': ObjectId(geneset_id)})
        if not geneset_result:
            bind_obj.set_error('can not find main_id: {} in sg_geneset'.format(geneset_id))
        task_id = geneset_result['task_id']
        geneset_type = geneset_result['type']
        geneset_name = geneset_result['name']
        genesets[geneset_name] = [geneset_type]
        names.append(geneset_name)
        collection = db['sg_geneset_detail']
        results = collection.find_one({'geneset_id': ObjectId(geneset_id)})
        geneset_names = set(results['seq_list'])
        genesets[geneset_name].append(geneset_names)
    # genesets.keys() == names
    return genesets, names, task_id, geneset_type

# class only
def export_cog_class(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    cog_path = os.path.join(dir_path, 'cog_class_table.xls')
    bind_obj.logger.debug('exporting {}'.format(cog_path))
    genesets, names, task_id, geneset_type = get_geneset_detail(data, bind_obj)
    cog_collection = db['sg_annotation_cog']
    cog_detail_collection = db['sg_annotation_cog_detail']
    cog_id = cog_collection.find_one({'task_id': task_id, 'type': bind_obj.sheet.option('type')})['main_id']
    one_record = cog_detail_collection.find_one({'cog_id': cog_id, 'anno_type': geneset_type})
    if not one_record:
        bind_obj.set_error('can not find cog_id: {} in sg_annotation_cog_detail'.format(cog_id))
    cog_results = cog_detail_collection.find({'cog_id': cog_id, 'anno_type': geneset_type})
    new_table_title = list()
    for tt in names:
        new_table_title.extend(['{}_COG'.format(tt), '{}_COG_list'.format(tt)])
    bind_obj.logger.debug(new_table_title)
    # write data in cog_class_table.xls
    with open(cog_path, 'wb') as w:
        w.write('Type\tFunctional Categoris\t{}\n'.format('\t'.join(new_table_title)))
        for cr in cog_results:
            if cr['cog_list']:
                cog_list = set(cr['cog_list'].split(';'))
            else:
                cog_list = list()
            write_line = dict()
            for gt in names:
                cog_count = list(cog_list & genesets[gt][1])
                if len(cog_count):
                    write_line[gt] = '{}\t{}'.format(len(cog_count),
                                                     ';'.join(cog_count))
            if len(write_line):
                w.write('{}\t{}\t'.format(cr['type'], cr['function_categories']))
                for tt in names:
                    if tt in write_line:
                        w.write('{}\t'.format(write_line[tt]))
                    else:
                        w.write('0\t\t')
                w.write('\n')
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, cog_path))
    return cog_path

# class only
def export_go_class(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    go_path = os.path.join(dir_path, 'go_class_table.xls')
    bind_obj.logger.info('exporting {}'.format(go_path))
    genesets, names, task_id, geneset_type = get_geneset_detail(data, bind_obj)
    go_collection = db['sg_annotation_go']
    go_level_collection = db['sg_annotation_go_detail']
    go_id = go_collection.find_one({'task_id': task_id, 'type': bind_obj.sheet.option('type')})['main_id']
    one_record = go_level_collection.find_one({'go_id': go_id, 'level': 2, 'anno_type': geneset_type})
    if not one_record:
        bind_obj.set_error('can not find go_id: {} in sg_annotation_go_detail'.format(go_id))
    new_table_title = list()
    for gt in genesets:
        new_table_title.append('{} num'.format(gt))
        new_table_title.append('{} percent'.format(gt))
        new_table_title.append('{} list'.format(gt))
    bind_obj.logger.debug(new_table_title)
    # write data in go_class_table.xls
    with open(go_path, 'wb') as w:
        w.write('Term type\tTerm\tGO\t{}\n'.format('\t'.join(new_table_title)))
        term_list = ['biological_process', 'cellular_component', 'molecular_function']
        for item in term_list:
            go_results = go_level_collection.find({'go_id': go_id, "level": 2, 'anno_type': geneset_type})
            for gr in go_results:
                if gr['goterm'] == item:
                    seq_list = set(gr['seq_list'].split(';'))
                    write_line = dict()
                    for gt in genesets:
                        total_gene_num = len(genesets[gt][1])
                        go_count = list(seq_list & genesets[gt][1])
                        if len(go_count):
                            write_line[gt] = '{}\t{}\t{}'.format(len(go_count),
                                                                 '{}({}/{})'.format(len(go_count)/total_gene_num,
                                                                                    len(go_count),
                                                                                    total_gene_num),
                                                                 ';'.join(go_count))
                    if len(write_line):
                        w.write('{}\t{}\t{}\t'.format(gr['goterm'], gr['goterm_2'], gr['goid_2']))
                        for tt in genesets:
                            if tt in write_line:
                                w.write('{}\t'.format(write_line[tt]))
                            else:
                                w.write('0\t0\t\t')
                        w.write('\n')
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, go_path))
    return go_path

# class and enrich
def export_kegg_table(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    kegg_table_path = os.path.join(dir_path, 'gene_kegg_table.xls')
    bind_obj.logger.debug('exporting {}'.format(kegg_table_path))
    genesets, names, task_id, geneset_type = get_geneset_detail(data, bind_obj)
    kegg_collection = db['sg_annotation_kegg']
    kegg_table_collection = db['sg_annotation_kegg_table']
    kegg_id = kegg_collection.find_one({'task_id': task_id, 'type': bind_obj.sheet.option('type')})['main_id']
    one_record = kegg_table_collection.find_one({'kegg_id': kegg_id, 'anno_type': geneset_type})
    if not one_record:
        bind_obj.set_error('can not find kegg_id: {} in sg_annotation_kegg_table'.format(ObjectId(kegg_id)))
    # write data in gene_kegg_table.xls
    with open(kegg_table_path, 'wb') as w:
        w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\n')
        kegg_results = kegg_table_collection.find({'kegg_id': kegg_id, 'anno_type': geneset_type})
        for result in kegg_results:
            if 'hyperlink' not in result:
                bind_obj.logger.debug('{} {} -> no hyperlink'.format(result['ko_id'], result['transcript_id']))
                result['hyperlink'] = 'None'
            w.write('{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths']))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, kegg_table_path))
    return kegg_table_path

# class only
def export_kegg_level(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    kegg_level_path = os.path.join(dir_path, 'gene_kegg_level_table.xls')
    bind_obj.logger.debug('exporting {}'.format(kegg_level_path))
    genesets, names, task_id, geneset_type = get_geneset_detail(data, bind_obj)
    kegg_collection = db['sg_annotation_kegg']
    kegg_level_collection = db['sg_annotation_kegg_level']
    kegg_id = kegg_collection.find_one({'task_id': task_id, 'type': bind_obj.sheet.option('type')})['main_id']
    one_record = kegg_level_collection.find_one({'kegg_id': kegg_id, 'anno_type': geneset_type})
    if not one_record:
        bind_obj.set_error('can not find kegg_id: {} in sg_annotation_kegg_level'.format(ObjectId(kegg_id)))
    my_result = kegg_collection.find({"task_id": task_id, "type": bind_obj.sheet.option("type")})
    # write data in gene_kegg_level_table.xls
    with open(kegg_level_path, 'wb') as w:
        w.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
        kegg_results = kegg_level_collection.find({'kegg_id': kegg_id, 'anno_type': geneset_type})
        for result in kegg_results:
            # graph_id (second column) and graph_png_id (ninth column) had beed deleted
            if 'hyperlink' not in result:
                bind_obj.logger.debug('{} -> no hyperlink'.format(result['pathway_id']))
                result['hyperlink'] = 'None'
            w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], # first
                                                                      '', # second
                                                                      result['number_of_seqs'], # third
                                                                      result['pathway_definition'], # fourth
                                                                      result['first_category'], # fifth
                                                                      result['anno_type'], # sixth
                                                                      result['hyperlink'], # seventh
                                                                      result['seq_list'], # eighth
                                                                      '', # ninth
                                                                      result['second_category'])) # tenth
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, kegg_level_path))
    return kegg_level_path

# class and enrich
def export_multi_gene_list(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    geneset_ids = data.split(',')
    multi_geneset_path = os.path.join(dir_path, 'multi_geneset_list.txt')
    bind_obj.logger.debug('exporting {}'.format(multi_geneset_path))
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    with open(multi_geneset_path, 'wb') as f:
        for n, gi in enumerate(geneset_ids):
            main_result = main_collection.find_one({'main_id': ObjectId(gi)})
            if not main_result:
                bind_obj.set_error('can not find main_id: {} in sg_geneset'.format(ObjectId(gi)))
            detail_result = collection.find_one({"geneset_id": ObjectId(gi)})
            if not detail_result:
                bind_obj.set_error('can not find geneset_id: {} in sg_geneset_detail'.format(ObjectId(gi)))
            f.write('{}\t{}\n'.format(main_result['name'], ','.join(detail_result['seq_list'])))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, multi_geneset_path))
    return multi_geneset_path

# class and enrich
def export_add_info(data, option_name,dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    task_id = data.split('\t')[0]
    anno_type = data.split('\t')[1]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug('exporting {}'.format(add_info))
    kegg_collection = db['sg_annotation_kegg']
    kegg_id = kegg_collection.find_one({'task_id': task_id, 'type': bind_obj.sheet.option('type')})["main_id"]
    kegg_level_collection = db['sg_annotation_kegg_level']
    kegg_results = kegg_level_collection.find({'kegg_id':kegg_id, 'anno_type': anno_type})
    with open(add_info, 'wb') as w:
        w.write('pathway\thyperlink\n')
        for result in kegg_results:
            w.write('{}\t{}\n'.format(result['pathway_id'], result['hyperlink']))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, add_info))
    return add_info

# enrich only
def export_gene_list(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    set_list_path = os.path.join(dir_path, 'geneset.list')
    bind_obj.logger.debug('exporting {}'.format(set_list_path))
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error('can not find main_id: {} in sg_geneset'.format(data))
    collection = db['sg_geneset_detail']
    results = collection.find_one({'geneset_id': ObjectId(data)})
    with open(set_list_path, 'wb') as f:
        gene_list = results['seq_list']
        for gene_id in gene_list:
            f.write('{}\n'.format(gene_id))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, set_list_path))
    return set_list_path

def export_geneset_list(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    set_list_path = os.path.join(dir_path, option_name + '.list')
    bind_obj.logger.debug('exporting {}'.format(set_list_path))
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error('can not find main_id: {} in sg_geneset'.format(data))
    collection = db['sg_geneset_detail']
    results = collection.find_one({'geneset_id': ObjectId(data)})
    with open(set_list_path, 'wb') as f:
        gene_list = results['seq_list']
        for gene_id in gene_list:
            f.write('{}\n'.format(gene_id))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, set_list_path))
    return set_list_path

# enrich only
def export_all_list(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    all_list_path = os.path.join(dir_path, 'all.list')
    bind_obj.logger.debug('exporting {}'.format(all_list_path))
    sg_geneset = db['sg_geneset']
    geneset_info = sg_geneset.find_one({'main_id': ObjectId(data)})
    if not geneset_info:
        bind_obj.set_error('can not find main_id: {} in sg_geneset'.format(data))
    task_id = geneset_info['task_id']
    project_sn = geneset_info['project_sn']
    all_seq = list()
    if bind_obj.sheet.option('anno_type') == 'go':
        sg_annotation_go = db['sg_annotation_go']
        main_table = sg_annotation_go.find_one({'task_id': task_id, 'project_sn': project_sn, 'type': bind_obj.sheet.option('type')})
        go_id = main_table['main_id']
        sg_annotation_go_detail = db['sg_annotation_go_detail']
        detail_tables = sg_annotation_go_detail.find({'go_id': ObjectId(go_id), 'anno_type': bind_obj.sheet.option('geneset_type')})
        for detail_table in detail_tables:
            all_seq.extend(detail_table['seq_list'].split(';'))
    elif bind_obj.sheet.option('anno_type') == 'kegg':
        sg_annotation_kegg = db['sg_annotation_kegg']
        main_table = sg_annotation_kegg.find_one({'task_id': task_id, 'project_sn': project_sn, 'type': bind_obj.sheet.option('type')})
        kegg_id = main_table['main_id']
        sg_annotation_kegg_level = db['sg_annotation_kegg_level']
        level_tables = sg_annotation_kegg_level.find({'kegg_id': ObjectId(kegg_id), 'anno_type': bind_obj.sheet.option('geneset_type')})
        for level_table in level_tables:
            all_seq.extend(level_table['seq_list'].split(';'))
    all_seq = set(all_seq)
    if len(all_seq):
        with open(all_list_path, 'wb') as f:
            for seq in all_seq:
                f.write('{}\n'.format(seq))
    else:
        bind_obj.set_error('no data in all_seq, anno_type: {}'.format(bind_obj.sheet.option('anno_type')))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, all_list_path))
    return all_list_path

# enrich only
def export_go_list(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    go_list_path = os.path.join(dir_path, 'go.list')
    bind_obj.logger.debug('exporting {}'.format(go_list_path))
    sg_geneset = db['sg_geneset']
    geneset_info = sg_geneset.find_one({'main_id': ObjectId(data)})
    if not geneset_info:
        bind_obj.set_error('can not find main_id: {} in sg_geneset'.format(data))
    task_id = geneset_info['task_id']
    project_sn = geneset_info['project_sn']
    sg_annotation_go = db['sg_annotation_go']
    main_table = sg_annotation_go.find_one({'task_id': task_id, 'project_sn': project_sn, 'type': bind_obj.sheet.option('type')})
    go_id = main_table['main_id']
    sg_annotation_go_list = db["sg_annotation_go_list"]
    list_tables = sg_annotation_go_list.find({'go_id': ObjectId(go_id)})
    with open(go_list_path, 'wb') as w:
        for table in list_tables:
            if table['anno_type'] == bind_obj.sheet.option('geneset_type'):
                w.write('{}\t{}\n'.format(table['gene_id'], table['gos_list']))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, go_list_path))
    return go_list_path

# ---------------------------------------- geneset related end -----------------------------------------

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
    group_schema = group_table.find_one({"main_id": ObjectId(data)})
    if not group_schema:
        bind_obj.set_error("无法根据传入的group_id:{}在sg_specimen_group_compare表里找到相应的记录".format(data))
    schema_name = re.sub("\s", "_", group_schema["group_name"])  # 将分组方案名的空格替换成下划线
    with open(file_path, "wb") as f:
        f.write("#sample\t" + "group" + "\n")
    sample_table_name = 'sg_specimen'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"main_id": ObjectId(sp_id)})
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
    file_path = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("正在导出计数矩阵:%s" % file_path)
    collection = db['sg_specimen_group_compare']
    result = collection.find_one({'main_id': ObjectId(data)})
    if not result:
        bind_obj.set_error("意外错误，control_id:{}在sg_specimen_group_compare中未找到！".format(ObjectId(data)))
    group_id = result['specimen_group_id']
    if group_id not in ['all', 'All', 'ALL']:
        """检查group_id的信息"""
        if isinstance(group_id, types.StringTypes):
            group_id = ObjectId(group_id)
        group_coll = db['sg_specimen_group']
        g_result = group_coll.find_one({'main_id': group_id})
        if not g_result:
            bind_obj.set_error("意外错误，control_file的group_id:{}在sg_specimen_group中未找到！".format(group_id))
    control_detail = json.loads(result['compare_names'])
    with open(file_path, 'wb') as w:
        w.write('#control\t{}\n'.format('group'))
        for i in control_detail:    #此处需要修改, 可能会有错误
            # w.write('{}\t{}\n'.format(i.keys()[0], i.values()[0]))
            control_other = i.split("|")
            w.write('{}\t{}\n'.format(control_other[0], control_other[1]))
    return file_path


def _get_gene_id(geneset, geneset_detail, _id, bind_obj=None):
    try:
        results = geneset_detail.find_one({"geneset_id":ObjectId(_id)})
        seq_id = results['seq_list']
    except Exception:
        bind_obj.set_error("{}在sg_geneset_detail表中没有找到!")
    try:
        my_result = geneset.find_one({"main_id":ObjectId(_id)})
        _name = my_result['name']
    except Exception:
        bind_obj.set_error("{}在sg_geneset表中没有找到!")
    return seq_id, _name

# 根据"add_info":geneset_info['task_id'] + "\t" + data.geneset_type，就是sg_geneset主表的task_id和页面传过来的是annotation还是annotation1
# 也就是是origin还是latest


# ------------------------------------------ exp related gdq -------------------------------------------

def export_exp_matrix(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    conn = db['sg_exp_detail']
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(seq_id=1, is_novel=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    exp_records = conn.find({'exp_id': ObjectId(data)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    ## modified by shicaiping @20190129， remove seq_id contain "." in plant novel miRNA
    remove_seq = exp_matrix[exp_matrix["seq_id"].str.contains(r"\.")]
    seq_id1 = list(remove_seq.seq_id)
    seq_id11 = list()
    for seq_id in seq_id1:
        if not "_" in seq_id.split(".")[1]:
            seq_id11.append(seq_id)
    seq_id2 = list(exp_matrix.seq_id)
    retain_id = list(set(seq_id2) ^ set(seq_id11))
    exp_matrix = exp_matrix[exp_matrix.seq_id.isin(retain_id)]

    # select specific seq_type for making correct exp_matrix
    if 'seq_type' in bind_obj.sheet.options():
        seq_type = bind_obj.sheet.option('seq_type')
        if seq_type.lower() == 'all':
            bind_obj.logger.debug('seq_type is {}, keep all miRNA'.format(seq_type))
        elif seq_type.lower() == 'known':
            bind_obj.logger.debug('seq_type is {}, keep known miRNA'.format(seq_type))
            exp_matrix = exp_matrix[exp_matrix['is_novel'] == False]
        elif seq_type.lower() == 'novel':
            bind_obj.logger.debug('seq_type is {}, keep novel miRNA'.format(seq_type))
            exp_matrix = exp_matrix[exp_matrix['is_novel'] == True]
    # end of filter
    exp_matrix.drop('is_novel', axis=1, inplace=True)
    exp_matrix = exp_matrix.set_index('seq_id')
    # make logarithmic transformation if provided use_log is True
    if 'use_log' in bind_obj.sheet.options() and bind_obj.sheet.option('use_log'):
        exp_matrix = np.log10(exp_matrix + 1)
    output = os.path.join(dir_path, '{}.xls'.format(option_name))
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, output))
    return output

def export_exp_matrix2(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    conn = db['sg_exp_detail']
    use_samples = bind_obj.sheet.option('use_samples').split("|")
    '''
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    '''
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in use_samples:
        target_cols[each] = 1
    exp_records = conn.find({'exp_id': ObjectId(data)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.loc[:, ["seq_id"] + use_samples]
    exp_matrix = exp_matrix.set_index('seq_id')
    output = os.path.join(dir_path, '{}.xls'.format(option_name))
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, output))
    return output

def export_group(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    group_dict = json.loads(data, object_pairs_hook=OrderedDict)
    output = os.path.join(dir_path, '{}.txt'.format(option_name))
    with open(output, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in group_dict:
            for each in group_dict[key]:
                f.write('{}\t{}\n'.format(each, key))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, output))
    return output

def export_compare(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    conn = db['sg_specimen_group_compare']
    result = conn.find_one({'main_id': ObjectId(data)})
    if not result:
        bind_obj.set_error('main_id: {} is not found in sg_specimen_group_compare'.format(ObjectId(data)))
    cmp_info = json.loads(result['compare_names'])
    output = os.path.join(dir_path, '{}.txt'.format(option_name))
    with open(output, 'w') as f:
        f.write('#ctrl\ttest\n')
        for each in cmp_info:
            f.write('{}\n'.format(each.replace('|', '\t')))
    bind_obj.logger.info('succeed in calling {} to return {}'.format(func_name, output))
    return output

# ------------------------------------------ exp related end -------------------------------------------

# ---------------------------------------- geneset related gdq -----------------------------------------

def export_geneset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数根据提供的geneset_id和exp_id和group_dict提取表达信息。
    当geneset_id为‘all’时，则导出所有的exp_id对应的详情表。
    当geneset_id为‘refall’时，则导出所有的exp_id对应的详情表时还需限制'is_new'字段为False(false)。
    该函数还顺便根据group_dict信息生成分组信息文件group_info.txt。
    该函数仅返回表达矩阵文件的路径信息
    """
    exp_id, geneset_id = data.split(";")
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
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    # get geneset
    if "all" not in geneset_id.lower():
        conn = db['sg_geneset_detail']
        geneset_records = conn.find_one({"geneset_id": ObjectId(geneset_id)})
        if not geneset_records:
            bind_obj.set_error('geneset not found by query: {}'.format(geneset_id))
        geneset = geneset_records['seq_list']
    else:
        geneset = list()
    # get all exp matrix
    conn = db['sg_exp_detail']
    if "refall" in geneset_id.lower():
        exp_records = conn.find({"exp_id": ObjectId(exp_id), 'is_new': False}, target_cols)
    else:
        exp_records = conn.find({"exp_id": ObjectId(exp_id)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    if exp_matrix.shape[0] == 0:
        bind_obj.set_error('No expression data find by query: {}'.format(exp_id))
    exp_matrix = exp_matrix.set_index('seq_id')
    if "all" not in geneset_id.lower():
        exp_matrix = exp_matrix.loc[geneset, :]
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_geneset_from_query(data, option_name, dir_path, bind_obj=None):
    func_name = sys._getframe().f_code.co_name
    chk_parm_func(func_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    main_collection_known = db["sg_known_mirna"]
    main_collection_novol = db["sg_novel_mirna"]
    collection_known = db["sg_known_mirna_detail"]
    collection_novol = db["sg_novel_mirna_detail"]
    collection = db["sg_annotation_query_detail"]
    main_collection = db["sg_annotation_query"]

    gene_type = bind_obj.sheet.option('gene_type')
    if gene_type == 'M':
        my_result_known = main_collection_known.find_one({'task_id': data})
        if not my_result_known:
            results_known = list()
            # 可能没有已知
            # bind_obj.set_error('task_id: {} is not found in sg_known_mirna'.format(data))
        else:
            if 'main_id' in my_result_known:
                query_id = my_result_known['main_id']
            else:
                query_id = my_result_known['_id']
            results_known = collection_known.find({'known_mirna_id': ObjectId(query_id)})

        my_result_novol = main_collection_novol.find_one({'task_id': data})

        if not my_result_novol:
            results_novol = list()
            # bind_obj.set_error('task_id: {} is not found in sg_novel_mirna'.format(data))
        else:
            if 'main_id' in my_result_novol:
                query_id = my_result_novol['main_id']
            else:
                query_id = my_result_novol['_id']
            results_novol = collection_novol.find({'novel_mirna_id': ObjectId(query_id)})

        output = os.path.join(dir_path, 'all_queries_id.txt')
        with open(output, 'w') as f:
            f.write('mirna\n')
            for result in results_known:
                mirna_id = ''
                if 'miRNA_name' in result:
                    mirna_id = result['miRNA_name']
                    f.write('{}\n'.format(mirna_id))
            if my_result_novol:
                for result in results_novol:
                    mirna_id = ''
                    if 'miRNA_id' in result:
                        mirna_id = result['miRNA_id']
                        f.write('{}\n'.format(mirna_id))
        return output

    elif gene_type == 'G':
        my_result = main_collection.find_one({'task_id': data})
        if not my_result:
            bind_obj.set_error('task_id: {} is not found in sg_annotation_query'.format(data))
        if 'main_id' in my_result:
            query_id = my_result['main_id']
        else:
            query_id = my_result['_id']
        results = collection.find({'query_id': ObjectId(query_id)})
        output = os.path.join(dir_path, 'all_queries_id.txt')
        with open(output, 'w') as f:
            f.write('gene\n')
            for result in results:
                gene_id = ''
                if 'gene_id' in result:
                    gene_id = result['gene_id']
                f.write('{}\n'.format(gene_id))
        return output

# ---------------------------------------- geneset related end -----------------------------------------

def export_wgcna_exp_matrix(data, option_name, dir_path, bind_obj=None):
    """ gdq
    该函数仅为wgcna module分析服务。
    该函数会返回wgcna_prepare得到的exp_matrix路径，
    同时还会通过查询sg_annotation_query表导出geneid和genename的关系文件:seq_id2gene_name.txt
    """
    wgcna_prepare = db['sg_wgcna_prepare']
    prepare_id = data.strip()
    prepare_main = wgcna_prepare.find_one({"_id": ObjectId(prepare_id)})
    exp_matrix = prepare_main['exp_matrix']
    os.link(exp_matrix, os.path.join(dir_path, "exp_matrix.txt"))
    target_seqs = pd.read_table(exp_matrix, header=0, index_col=0).index
    task_id = prepare_main['task_id']
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "latest"})
    if not annot_main:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_annotation_query_detail']
    query_dict = dict(query_id=annot_main_id,)
    result_dict = dict( _id=0, gene_name=1, gene_id=1, transcript_id=1)
    result = annot_detail.find(query_dict, result_dict)
    gene2name = pd.DataFrame(list(result))
    exp_level = bind_obj.sheet.option('exp_level')
    if exp_level[0].upper() == 'T':
        gene2name.set_index('transcript_id', inplace=True)
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
    module_id = data.strip()
    # export eigengenes
    eigengenes = db['sg_wgcna_module_eigengenes_detail']
    eigengenes_found = eigengenes.find({"module_id": ObjectId(module_id)}, {"_id": 0, "module_id":0})
    eigengenes_pd = pd.DataFrame(list(eigengenes_found))
    eigengenes_pd.set_index("module", inplace=True)
    eigengenes_path = os.path.join(dir_path, "module_eigengenes.xls")
    eigengenes_pd.to_csv(eigengenes_path, sep='\t', header=True, index=True)
    # export exp matrix
    prepare_id = db['sg_wgcna_module'].find_one({"main_id": ObjectId(module_id)})['wgcna_prepare_id']
    exp_matrix = db['sg_wgcna_prepare'].find_one({"main_id": ObjectId(prepare_id)})['exp_matrix']
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


def get_after_qc_bam_path(data, option_name, dir_path, bind_obj=None):
    db = Config().get_mongo_client(mtype="smallrna")[Config().get_mongo_dbname("smallrna")]
    conn = db['sg_specimen']
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
    conn = db['sg_diff_detail']
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
    conn = db['sg_geneset_go_enrich_detail']
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
    conn = db['sg_geneset_kegg_enrich_detail']
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
    collection = db['sg_geneset_detail']
    main_collection = db['sg_geneset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("意外错误，geneset_id:{}在sg_geneset中未找到！".format(ObjectId(data)))
    results = collection.find_one({"geneset_id": ObjectId(data)})["seq_list"]
    with open(gene_list_path, "wb") as f:
        f.write("accession_id" + "\n")
        for result in results:
            f.write(result + "\n")
    bind_obj.logger.debug("蛋白集导出成功！")
    return gene_list_path

# added for transcription factor analysis  --------------------------------------------
def get_all_pep_seq(data, option_name, dir_path, bind_obj=None):
    pep_db_path = data.strip()
    result_path = os.path.join(dir_path, "all_pep.fa")
    import sqlite3
    conn = sqlite3.connect(pep_db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT transcript_id,pep_seq FROM trans_annot")
    with open(result_path, 'w') as fw:
        # pep_id 就是转录本的id
        for pep_id, pep_seq in cursor.fetchall():
            fw.write('>{}\n{}\n'.format(pep_id, pep_seq))
    return result_path


def get_gene_detail(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip()
    annot_table = db['sg_annotation_query']
    annot_main = annot_table.find_one({"task_id": task_id, "type": "latest", 'status': 'end'})
    if not annot_main:
        annot_main = annot_table.find_one({"task_id": task_id, "type": "origin"})
    if not annot_main:
        bind_obj.set_error("Not Found in sg_annotation_query by query {}".format(task_id))
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    annot_detail = db['sg_annotation_query_detail']
    query_dict = dict(query_id=annot_main_id, )
    result_dict = dict(_id=0, gene_name=1, gene_id=1, description=1, transcript_id=1)
    result = annot_detail.find(query_dict, result_dict)
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

def get_gene_link(data, option_name, dir_path, bind_obj=None):
    task_id = data.strip()
    annot_table = db['sg_genes']
    annot_main = annot_table.find_one({"task_id": task_id})
    if not annot_main:
        gene_lnc = os.path.join(dir_path, "gene_lnc.xls")
        with open(gene_lnc, 'w') as f:
            f.write("gene_id\tgene_ensembl")
        return gene_lnc
        # bind_obj.set_error("Not Found in sg_annotation_query by query {}".format(task_id))
    if "main_id" not in annot_main:
        annot_main_id = annot_main['_id']
    else:
        annot_main_id = annot_main['main_id']
    lnc_detail = db['sg_genes_detail']
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




# added by gdq for tfbs predict
def export_geneid2tfid_file(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_tf_predict_detail']
    result = conn.find({"tf_predict_id": ObjectId(data.strip())}, {"_id": 0, 'gene_id': 1, "blast_hit": 1})
    result_pd = pd.DataFrame(list(result))
    result_pd.set_index("gene_id", inplace=True)
    result_file = os.path.join(dir_path, "geneid2tfid.txt")
    result_pd.to_csv(result_file, sep='\t', header=True, index=True)
    return result_file


def export_predict_result(data, option_name, dir_path, bind_obj=None):
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

# ---------------------------------------------------------------------------------------------

# ---------------------------------------- check parameter qjc -----------------------------------------

def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} - {}'.format(k, v))
    else:
        pass

def chk_parm_dcrt(func):
    '''
    this decorator was deprecated, because the use of decorator resulted in returning None of to_file func
    '''
    def wrapper(*args, **kwargs):
        if 'bind_obj' in kwargs:
            kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func.__name__))
            for i, e in enumerate(args):
                kwargs['bind_obj'].logger.debug('index_{} - {}'.format(i, e))
            for k, v in kwargs.iteritems():
                kwargs['bind_obj'].logger.debug('{} - {}'.format(k, v))
        else:
            for i, e in enumerate(args):
                if hasattr(e, 'id'):
                    wf_loc = i
            args[wf_loc].logger.info('check to_file parameters in {}'.format(func.__name__))
            for i, e in enumerate(args):
                args[wf_loc].logger.debug('index_{} - {}'.format(i, e))
            for k, v in kwargs.iteritems():
                args[wf_loc].logger.debug('{} - {}'.format(k, v))
        func(*args, **kwargs)
    return wrapper

# ---------------------------------------- check parameter end -----------------------------------------

#新增批次效应的to_file函数(export_exp_batch_matrix, export_group_detail, export_exp_other_matrix)

def export_exp_batch_matrix(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    exp_id, sample_str = data.split(';')
    sample_list = sample_str.split(',')
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in sample_list:
        target_cols[each] = 1
    collection = db['sg_exp_detail']
    print target_cols
    cursor = collection.find({'exp_id': ObjectId(exp_id)}, target_cols)
    output = os.path.join(dir_path, 'count.txt')
    df = pd.DataFrame(list(cursor))
    exp_matrix = df.set_index('seq_id')
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    # df.to_csv(output, sep='\t', index=False)
    return output

def export_group_detail(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    collection = db['sg_specimen_group']
    record = collection.find_one({'task_id': data})
    group = record['category_names']
    sample = record['specimen_names']
    group_dict = OrderedDict()
    for i in range(len(group)):
        group_dict[group[i]] = sample[i]
    output = os.path.join(dir_path, 'group_table.txt')
    with open(output, 'w') as f:
        f.write('#sample\tgroup\n')
        for each in group_dict:
            for key in group_dict[each]:
                f.write('{}\t{}\n'.format(key, each))
    return output

def export_exp_other_matrix(data, option_name, dir_path, bind_obj):
    chk_parm_func(
        sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    )
    exp_id = data
    collection = db['sg_exp_detail']
    # if exp_level == 'T':
    #     other_cols = OrderedDict(seq_id=1, _id=0, is_new=1, gene_id=1)
    # else:
    #     other_cols = OrderedDict(seq_id=1, _id=0, is_new=1)
    other_cols = OrderedDict(seq_id=1, _id=0, is_novel=1)
    other = collection.find({'exp_id': ObjectId(exp_id)}, other_cols)
    output_other = os.path.join(dir_path, 'other.txt')
    df_other = pd.DataFrame(list(other))
    other_matrix = df_other.set_index('seq_id')
    other_matrix.to_csv(output_other, sep='\t', header=True, index=True)
    # df.to_csv(output, sep='\t', index=False)
    return output_other
