# -*- coding: utf-8 -*-
# __author__ = 'HongDong'
# modified = 20180919

import os
import re
import json
from bson import SON
from .core.base import Base
import datetime
from bson.objectid import ObjectId


class Dna(Base):
    """
    该脚本用于dna项目的接口需要的各种数据库操作，包含查找与更新，可以根据各自项目在后面继续完善
    modified by hongdong 20180409
    """
    def __init__(self, db_types):
        super(Dna, self).__init__()
        self._project_type = db_types

    def find_one(self, collection, query_dic):
        collection = self.db[collection]
        return collection.find_one(query_dic)

    def find_many(self, collection, query_dic):
        collection = self.db[collection]
        return collection.find(query_dic)

    def insert_one(self, collection, data_dic):
        collection = self.db[collection]
        return collection.insert_one(data_dic).inserted_id

    def insert_many(self, collection, data_list):
        collection = self.db[collection]
        collection.insert_many(data_list)

    def insert_main_table(self, collection, data):
        """
        用于插入主表
        :param collection: 集合的名称
        :param data: 插入内容
        :return:
        """
        try:
            id_ = self.db[collection].insert_one(SON(data)).inserted_id
        except Exception as e:
            raise Exception("插入主表{}失败！{}".format(collection, e))
        else:
            print "插入主表{}成功！".format(collection)
        return id_

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
            file_path_ = region.rstrip(":") + ':' + file_path
        elif file_path.startswith("/mnt") or re.match(".*://.*", file_path):
            file_path_ = file_path
        else:
            raise Exception("存入的文件{}路径格式不正确！".format(file_path))
        return file_path_

    def set_main_table_name(self, analysis_name, chongmingming_result, task_id=None, collection=None):
        """
        因为我们接口要添加自定义文件结果的功能，所以我们对结果目录main_table_name名字进行重新定义一下
        这里要注意一下文件夹名字的话是否会有重复，理论上李小杰那边进行了判断，我们这边暂时不进行判断
        :param analysis_name:  PopTree
        :param chongmingming_result:
        :param task_id:
        :param collection:
        :return:
        """
        if task_id and collection:
            result = self.db[collection].find({"task_id": task_id})
            for m in result:
                if m['name'] == chongmingming_result:
                    raise Exception("文件结果名**{}**已经存在，请重新输入其它文件夹名称！".format(chongmingming_result))
        if chongmingming_result:
            main_table_name = chongmingming_result + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        else:
            main_table_name = analysis_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        return main_table_name
