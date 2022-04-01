# !/usr/bin/python
# -*- coding: utf-8 -*-


import types
import os
import json
import unittest
import datetime
import sqlite3
import re
from collections import OrderedDict, defaultdict
from bson.son import SON
from bson.objectid import ObjectId
from mbio.api.database.prok_rna.api_base import ApiBase
from biocluster.api.database.base import Base, report_check
import glob
import copy


class TfPrediction(ApiBase):
    def __init__(self, bind_object):
        super(TfPrediction, self).__init__(bind_object)
        self._project_type = 'prok_rna'

    def add_tf_detail(self, tf_id, result):
        if not os.path.exists(result):
            raise Exception('{}所指定的路径不存在，请检查！'.format(result))
        data_list = list()
        if os.path.exists(result):
            with open(result, 'r') as f:
                f.readline()
                for line in f:
                    line = line.strip('\n').split('\t')
                    data = [
                        ('tf_id', ObjectId(tf_id)),
                        ('gene_id', line[0]),
                        ('class', line[1]),
                        ('type', line[2]),
                        ('p2tf_des', line[3]),
                        ('domain_id', line[4]),
                        ('domain_arch', line[5]),
                        ('identity', float(line[6])),
                        ('coverage', round(float(line[7])*100, 1)),
                        ('evalue', float(line[8])),
                        ('score', float(line[9])),
                    ]
                    data = SON(data)
                    data_list.append(data)
            try:
                self.create_db_table('sg_tf_predict_detail', data_list)
                self.db['sg_tf_predict'].update({"_id": ObjectId(tf_id)}, {
                    "$set": {"main_id": ObjectId(tf_id)}})
            except:
                self.bind_object.logger.error("导入转录因子注释信息：%s出错!" % (result))
            else:
                self.bind_object.logger.info("导入转录因子注释信息：%s 成功!" % result)