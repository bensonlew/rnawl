# -*- coding: utf-8 -*-
# __author__ = 'sanger'

from biocluster.config import Config
import os
import json
from bson.objectid import ObjectId
from collections import OrderedDict
import sys
import pandas as pd


project_type = 'ref_rna_v2'
db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

def checkwargs(**kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('start checking to_file arguments')
        for k, v in kwargs.items():
            kwargs['bind_obj'].logger.debug('{} = {}'.format(k, v))

def export_bam_list(data, option_name, dir_path, bind_obj=None):
    checkwargs(data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    results = db['sg_specimen'].find({'task_id': data, 'about_qc': 'after'})
    output = os.path.join(dir_path, 'bam.list')
    open(output, 'w').writelines(sorted(['{}\n'.format(i['bam_path']) for i in results]))
    return output

def export_rmats_group_table(data, option_name, dir_path, bind_obj):
    checkwargs(data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    lines = ['#sample\tgroup\n']
    for group, samples in json.loads(data).items():
        lines.extend(['{}\t{}\n'.format(sample, group) for sample in sorted(samples)])
    group_table = os.path.join(dir_path, 'group.txt')
    open(group_table, 'w').writelines(lines)
    return group_table

def export_rmats_control_table(data, option_name, dir_path, bind_obj):
    checkwargs(data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    control_table = os.path.join(dir_path, 'control.txt')
    open(control_table, 'w').writelines(['#control\tother\n', '{}\t{}\n'.format(*data.split('|'))])
    return control_table

def export_rmats_root(data, option_name, dir_path, bind_obj):
    checkwargs(data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    return '{}/'.format(db['sg_splicing_rmats'].find_one({'main_id': ObjectId(data)})['result_dir'])

def export_rmats_detail_path2base(data, option_name, dir_path, bind_obj):
    checkwargs(data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    lines = list()
    if 'task_id' in bind_obj.sheet.data:
        task_id = bind_obj.sheet.data['task_id']
    else:
        task_id = "_".join(bind_obj.sheet.data['id'].split("_")[0:-2])
    for splicing_id in data.split(','):
        rmats_info = db['sg_splicing_rmats'].find_one(
            {'main_id': ObjectId(splicing_id), 'task_id': task_id}
        )
        lines.append('{}\t{}\n'.format(
            os.path.join(rmats_info['rmats_output'], 'all_events_detail_big_table.txt'),
            '{}.txt'.format(rmats_info['main_id'])
        ))
    output = os.path.join(dir_path, 'table.list')
    open(output, 'w').writelines(lines)
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


def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} - {}'.format(k, v))
    else:
        pass

def export_fusion_matrix(data, option_name, dir_path, bind_obj=None):

    fusion_id = bind_obj.sheet.option('fusion_id')
    target_cols = OrderedDict(splicetype=1, leftbreakpoint=1, fusionname=1, spanningfragcount=1, rightbreakpoint=1,fusion_unique_id=1,
                              junctionreadcount=1, rightgene=1, leftgene=1, sample=1, _id=0)
    # target_cols = OrderedDict(SpliceType=1, LeftBreakpoint=1, FusionName=1, SpanningFragCount=1, RightBreakpoint=1, JunctionReadCount=1, RightGene=1,LeftGene=1, sample=1,_id=0)
    bind_obj.logger.debug("导出fusion分析表")
    # get geneset
    conn = db['sg_gene_fusion_detail']

    fusion_records = conn.find({"gene_fusion_id": ObjectId(fusion_id)}, target_cols)
    fusion_matrix = pd.DataFrame(list(fusion_records))
    output = os.path.join(dir_path, option_name)
    fusion_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export fusion_detail matrix')
    return output
