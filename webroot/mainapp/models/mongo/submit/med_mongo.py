# -*- coding: utf-8 -*-
# __author__ = 'yuguo'
from mainapp.models.mongo.core.base import Base
from biocluster.config import Config
from mainapp.libs.param_pack import param_pack
from bson import ObjectId
from bson import SON
import datetime


class MedMongo(Base):
    """
    医学数据拆分与亲子鉴定、产前筛查等流程相关数据库操作。(V2.version)
    modified by hongdong 20171121 亲子与产筛是不同数据库，添加不同库的连接
    """
    def __init__(self, db_types):
        super(MedMongo, self).__init__()
        self._project_type = db_types
        # if db_types == 'pt_v2':
        #     self._project_type = 'pt_v2'
        # else:
        #     self._project_type = 'nipt_v2'

    def get_one_bymore(self, collection, query_dic):
        collection = self.db[collection]
        return collection.find_one(query_dic)

    def get_many_bymore(self, collection, query_dic):
        collection = self.db[collection]
        return collection.find(query_dic)

    def get_one(self, collection, query, value):
        collection = self.db[collection]
        return collection.find_one({query: value})

    def get_many(self, collection, query, value):
        collection = self.db[collection]
        return collection.find({query: value})

    def add_one(self, collection, data_dic):
        collection = self.db[collection]
        return collection.insert_one(data_dic).inserted_id

    def add_many(self, collection, data_list):
        collection = self.db[collection]
        return collection.insert_many(data_list).inserted_ids

    def add_main_id(self, collection):
        return self.db[collection].insert_one({}).inserted_id

    def update_table(self, collection, query, value, dict):
        """
        用于更新一个表的一个字段或者多个字段
        :param collection:
        :param query:
        :param value:
        :param dict:
        :return:
        """
        try:
            self.db[collection].update({query: value}, {'$set': dict}, upsert=True, multi=True)
        except Exception as e:
            raise Exception("更新{}失败！{}".format(collection, e))
