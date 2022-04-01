# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from mainapp.config.db import get_mongo_client
from bson.objectid import ObjectId
import types
from biocluster.config import Config
from mainapp.models.workflow import Workflow
from meta import Meta


class Denovo(Meta):
    def __init__(self, db=None):
        self.client = get_mongo_client()
        if not db:
            self.db_name = Config().MONGODB + '_rna'
        else:
            self.db_name = db
        self.db = self.client[self.db_name]

    def get_main_info(self, main_id, collection_name):
        if isinstance(main_id, types.StringTypes):
            main_id = ObjectId(main_id)
        elif isinstance(main_id, ObjectId):
            main_id = main_id
        else:
            raise Exception("输入main_id参数必须为字符串或者ObjectId类型!")
        collection = self.db[collection_name]
        express_info = collection.find_one({'_id': main_id})
        return express_info
