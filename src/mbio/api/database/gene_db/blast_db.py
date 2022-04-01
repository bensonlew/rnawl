# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
import os
import re
import sys
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
import gridfs
import json
import unittest
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
from api_base import ApiBase


class BlastDb(ApiBase):
    def __init__(self, bind_object):
        super(BlastDb, self).__init__(bind_object)


    def add_blastpath_db(self, blast_path):
        '''
        添加数据库物名称信息
        '''
        collection = self.db["blast_db"]
        data_list = []
        with open(blast_path, 'r') as f:
            for line in f.readlines():
                cols = line.split("\t")
                data = {
                    "name": cols[0],
                    "version": cols[1],
                    "path": cols[2]
                }
                data_list.append(SON(data))
        if data_list:
            try:
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入blast信息失败" )




if __name__ == '__main__':
    BlastDb(None).add_blastpath_db(sys.argv[1])
