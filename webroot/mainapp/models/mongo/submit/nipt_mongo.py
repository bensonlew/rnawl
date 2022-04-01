# -*- coding: utf-8 -*-
# __author__ = 'hongdong.xuan'
from bson import SON
from bson import ObjectId
from biocluster.config import Config
from types import StringTypes
from mainapp.models.mongo.core.base import Base


class NiptMongo(Base):
    '''
    获取流程的样本名
    '''

    def __init__(self):
        # self.mongo_client = Config().mongo_client
        # self.database = self.mongo_client[Config().MONGODB + '_nipt']
        super(NiptMongo, self).__init__()
        self._project_type = "nipt"
        # self.client = Config().get_mongo_client(mtype=self._project_type)
        # self.database = self.client[Config().get_mongo_dbname(self._project_type)]

    def get_sample_info(self, _id):
        if not isinstance(_id, ObjectId):
            if isinstance(_id, StringTypes):
                _id = ObjectId(_id)
            else:
                raise Exception('_id必须为ObjectId对象或其对应的字符串!')
        # collection = self.database['sg_main']
        collection = self.db['sg_main']
        task_info = collection.find_one({"_id": _id})
        return task_info

    def insert_main_table(self, collection, data):
        # return self.database[collection].insert_one(SON(data)).inserted_id
        return self.db[collection].insert_one(SON(data)).inserted_id

    def insert_none_table(self, collection):
        # return self.database[collection].insert_one({}).inserted_id
        return self.db[collection].insert_one({}).inserted_id
