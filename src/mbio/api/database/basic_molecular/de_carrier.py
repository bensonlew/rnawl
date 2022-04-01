# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'

from collections import OrderedDict
from bson import SON
from bson.objectid import ObjectId
from biocluster.config import Config
from types import StringTypes
from io import StringIO
import types
import re
import datetime
import os
import json
import time
# from api_base import ApiBase
from biocluster.api.database.base import Base, report_check

class DeCarrier(Base):
    def __init__(self, bind_object):
        super(DeCarrier, self).__init__(bind_object)
        self._project_type = 'basic_molecular'

    def update_decarrier_record(self, main_id, s3_upload_path):
        """
        更新主表字段
        """
        update_dict = {
            "result_path": s3_upload_path,
        }
        self.db['sg_decarrier'].update({"_id":self.check_objectid(main_id)},{'$set':update_dict}, upsert =True, multi= True)


    
    
    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        :param id_:
        :return:
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_