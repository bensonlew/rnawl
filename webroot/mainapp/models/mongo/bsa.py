# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified = 20180223

import os
import json
import re
from bson import SON
from .core.base import Base
from bson.objectid import ObjectId


class Bsa(Base):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(Bsa, self).__init__(self._bind_object)
        self._project_type = "bsa"

    def insert_main_table(self, collection, data):
        """
        用于插入主表
        :param collection: 集合的名称
        :param data: 插入内容
        :return:
        """
        try:
            return self.db[collection].insert_one(SON(data)).inserted_id
        except Exception as e:
            raise Exception("插入主表{}失败！{}".format(collection, e))

    def update_db_record(self, collection, query_dict, update_dict, upsert=False, multi=False):
        """
        用于更新表格的字段，暂时只是简单的更新一个表的对应字段，后面有需求再进行完善
        :param collection: 集合的名称
        :param query_dict: 查询的字段，可以有多个字段进行联合查询, 是一个字典
        :param update_dict: 表格中要更新进去的字段
        :param upsert: 表中没有查找到对应条件的记录，确认是否要更新 默认为True
        :param multi: 表中按照对应条件查找到多条记录，是否要全部更新 默认为True
        :return:
        """
        try:
            self.db[collection].update(query_dict, {'$set': update_dict}, upsert=upsert, multi=multi)
        except Exception as e:
            raise Exception("更新表格{}失败！{}".format(collection, e))
        else:
            print "更新表格{}成功！".format(collection)

    def find_one_record(self, collection, query_dic):
        """
        用于mongo查询一条记录，查询字段是一个字典,支持一个字段的查询，也支持多个字段的查询
        example: {"batch_id": ObjectID("5a782577a4e1af477e1ac081"), "type": "pt"}
        :param collection:
        :param query_dic:
        :return:
        """
        collection = self.db[collection]
        return collection.find_one(query_dic)

    def find_records(self, collection, query_dic):
        """
        用于mongo查询多条记录，查询字段是一个字典,支持一个字段的查询，也支持多个字段的查询
        example: {"batch_id": ObjectID("5a782577a4e1af477e1ac081"), "type": "pt"}
        :param collection:
        :param query_dic:
        :return:
        """
        collection = self.db[collection]
        return collection.find(query_dic)

    def set_file_path(self, task_id, file_path, client):
        """
        该函数是用于封装一下所有的文件的路径的重构，包含了/mnt磁盘路径，以及rere开始的路径，以及//开始的S3路径
        add by hongdong@20180919
        :param task_id:
        :param file_path:
        :param client:
        :return:
        """
        if file_path.startswith("rerewrweset"):
            base_path = '/mnt/ilustre/data/' if str(client) == 'client01' else "/mnt/ilustre/tsanger-data/"
            # file_path_ = os.path.join(base_path, file_path)
            file_path_ = base_path + file_path
        elif file_path.startswith("//"):
            result = self.db["sg_task"].find_one({"task_id": task_id})
            # noinspection PyBroadException
            try:
                region = result['region']  # 's3:'
            except:
                raise Exception("sg_task中没有找到region字段！")
            # file_path_ = os.path.join(region, file_path)  # 不知道为啥不行
            file_path_ = region.rstrip(':') + ":" + file_path
        elif file_path.startswith("/mnt") or re.match(".*://.*", file_path):
            file_path_ = file_path
        else:
            raise Exception("存入的文件{}路径格式不正确！".format(file_path))
        return file_path_
