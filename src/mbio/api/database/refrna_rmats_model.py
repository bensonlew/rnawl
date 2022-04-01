# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/17 11:36

import re, os, Bio, argparse, sys, fileinput, urllib2

from pymongo import MongoClient
from bson.objectid import ObjectId
import types
from types import StringTypes
import re, subprocess
import json, time
import pandas as pd
import numpy as np
import datetime, os
from bson.son import SON
from collections import Counter
import glob
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import gridfs
from mainapp.models.mongo.ref_rna import RefRna

'''
rmats model 导表函数
'''


class RefrnaRmatsModel(Base):
    def __init__(self, bind_object):
        super(RefrnaRmatsModel, self).__init__(bind_object)
        self._project_type = 'ref_rna'
        #self._db_name = Config().MONGODB + '_ref_rna'
    
    @report_check
    def add_rmats_model(self, params=None, splicing_id=None, name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        splicing_id = splicing_id
        data = {'splicing_id': splicing_id,
                'name': name if name else 'rmats_model_graph_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
                'desc': '可变剪接rmats事件模型图主表',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'params':
                    json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params,
                'status': 'end',
                'task_id': task_id,
                'project_sn':project_sn
                }
        self.bind_object.logger.info('要插入的数据为：%s' % data)
        collection_obj = self.db['sg_splicing_rmats_model']
        try:
            rmats_model_id = collection_obj.insert_one(data).inserted_id
            self.bind_object.logger.info('导入rmats model主表成功')
        except Exception as e:
            raise Exception('导入rmats model主表失败:%s' % e)
        return rmats_model_id
    
    def add_sg_fs(self, file_path, rmats_model_id,file_type):
        data = open(file_path, "r").read()
        fs = gridfs.GridFS(database=self.db, collection='fs')
        fs_id = fs.put(data, filename=os.path.basename(file_path))  # splicing_id
        # self.db['fs.files'].update_one({'_id': fs_id}, {'$set': {'splicing_id': splicing_id}})
        # self.db['fs.files'].update_one({'_id': fs_id}, {'$set': {'rmats_model_id': rmats_model_id}})
        if file_type == 'pdf':
            self.db['sg_splicing_rmats_model'].update_one({'_id': rmats_model_id}, {'$set': {'graph_id': fs_id}})
        if file_type == 'png':
            self.db['sg_splicing_rmats_model'].update_one({'_id': rmats_model_id}, {'$set': {'png_id': fs_id}})

