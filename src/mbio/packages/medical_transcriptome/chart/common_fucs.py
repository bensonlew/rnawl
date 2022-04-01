# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import glob
import os
import re
import shutil
import types
import numpy as np
import pandas as pd
from biocluster.config import Config
from bson.objectid import ObjectId



# self.database = Config().get_mongo_client(mtype='medical_transcriptome')[Config().get_mongo_dbname('medical_transcriptome')]

class CommonFucs(object):
    def __init__(self,project_type="medical_transcriptome"):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

    def get_main_info(self, main_id, collection_name, task_id):
        """
        根据主表id获得整个主表/记录的信息
        :param main_id:
        :param collection_name:
        :return: 整个主表信息
        """
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db[collection_name]
        main_info = collection.find_one({'main_id': main_id, 'task_id': task_id})
        return main_info

    def get_rmats_stats_params_info(self, main_id):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db['sg_splicing_rmats_stats']
        result = collection.find_one({'main_id': main_id})
        return result