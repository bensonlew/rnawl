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

def checkwargs(**kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('start checking to_file arguments')
        for k, v in kwargs.items():
            kwargs['bind_obj'].logger.debug('{} = {}'.format(k, v))

def export_rmats_root(data, option_name, dir_path, bind_obj):
    checkwargs(data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    return '{}/'.format(db['splicing_rmats'].find_one({'main_id': ObjectId(data)})['result_dir'])

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

def export_bam_list(data, option_name, dir_path, bind_obj=None):
    checkwargs(data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    results = db['specimen'].find({'task_id': data, 'about_qc': 'after'})
    output = os.path.join(dir_path, 'bam.list')
    open(output, 'w').writelines(sorted(['{}\n'.format(i['bam_path']) for i in results]))
    return output

def export_rmats_detail_path2base(data, option_name, dir_path, bind_obj):
    checkwargs(data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj)
    lines = list()
    if 'task_id' in bind_obj.sheet.data:
        task_id = bind_obj.sheet.data['task_id']
    else:
        task_id = "_".join(bind_obj.sheet.data['id'].split("_")[0:-2])
    for splicing_id in data.split(','):
        rmats_info = db['splicing_rmats'].find_one(
            {'main_id': ObjectId(splicing_id), 'task_id': task_id}
        )
        lines.append('{}\t{}\n'.format(
            os.path.join(rmats_info['rmats_output'], 'all_events_detail_big_table.txt'),
            '{}.txt'.format(rmats_info['main_id'])
        ))
    output = os.path.join(dir_path, 'table.list')
    open(output, 'w').writelines(lines)
    return output