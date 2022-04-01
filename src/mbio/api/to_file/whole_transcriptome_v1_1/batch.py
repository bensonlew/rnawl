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

def chk_parm_func(func_name, **kwargs):
    if 'bind_obj' in kwargs and hasattr(kwargs['bind_obj'], 'id'):
        kwargs['bind_obj'].logger.info('check to_file parameters in {}'.format(func_name))
        for k, v in kwargs.iteritems():
            kwargs['bind_obj'].logger.debug('{} - {}'.format(k, v))
    else:
        pass


def export_exp_matrix_all(data, option_name, dir_path, bind_obj):
    # chk_parm_func(
    #     sys._getframe().f_code.co_name, data=data, option_name=option_name, dir_path=dir_path, bind_obj=bind_obj
    # )
    exp_id, sample_str, exp_level, is_rmbe = data.split(';')
    sample_batch_list = ['{}_batch'.format(x) for x in sample_str.split(',')]
    target_cols = OrderedDict(_id=0)
    # for each in sample_batch_list:
    #     target_cols[each] = 0
    collection = db['exp_detail']
    print target_cols
    cursor = collection.find({'exp_id': ObjectId(exp_id), }, target_cols)
    output = os.path.join(dir_path, 'count_all.txt')
    df = pd.DataFrame(list(cursor))
    # if exp_level == 'T':
    #     df.rename(columns={'transcript_id': 'seq_id'}, inplace=True)
    # if exp_level == 'G':
    #     df.rename(columns={'gene_id': 'seq_id'}, inplace=True)
    if exp_level == 'G':
        exp_matrix = df.set_index('gene_id')
    else:
        exp_matrix = df.set_index('transcript_id')
    exp_matrix.to_csv(output, sep='\t', header=True,index=True)
    # df.to_csv(output, sep='\t', index=False)
    return output