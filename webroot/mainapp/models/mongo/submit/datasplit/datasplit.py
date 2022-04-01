# -*- coding: utf-8 -*-
# __author__ = "wangzhaoyue"

import json
import datetime
from bson import SON
from bson.objectid import ObjectId
from mainapp.models.mongo.core.base import Base


class Datasplit(Base):
    def __init__(self, db_types="datasplit"):
        super(Datasplit, self).__init__()
        self._project_type = db_types

    def get_bcl2fastq_info(self, status_id):
        collection = self.db["seq_status"]
        print collection
        result = collection.find_one({"_id": ObjectId(status_id)})
        print result
        board_id = result["seq_board_id"]
        collection2 = self.db["seq_board"]
        board_result = collection2.find_one({"_id": board_id})
        board_data_path = board_result["data_path"]
        seq_number = board_result["seq_number"]
        try:
            params = json.loads(result["params"])
            board_seq_model = params["library_split"]["seq_model"]
            if "barcode_mismatch" in params["library_split"].keys():
                barcode_mismatch = int(params["library_split"]["barcode_mismatch"])
            else:
                barcode_mismatch = 0
        except:
            board_seq_model = board_result["seq_model"]
        return board_data_path, board_seq_model, seq_number, barcode_mismatch

    def coll_find_one(self, collection, query_dict):
        """
        查找表中符合条件的一条记录
        collection: 集合的名称
        query_dict: 查找条件，是一个字典
        """
        result = self.db[collection].find_one(query_dict)
        return result

    def coll_find(self, collection, query_dict):
        """
        查找表中符合条件的所有记录
        collection: 集合的名称
        query_dict: 查找条件，是一个字典
        """
        result = self.db[collection].find(query_dict)
        return result

    def update_db_record(self, collection, query_dict, update_dict, is_show_log="true", upsert=True, multi=True):
        """
        用于更新表格的字段
        collection: 集合的名称
        query_dict: 查找条件，是一个字典
        update_dict: 表格中要更新进去的字段，是一个字典
        is_show_log: 用于设定是否要显示出更新成功的日志信息, 为false的时候就不显示成功的日志信息
        upsert: 表中没有查找到对应条件的记录，确认是否要更新 默认为True
        multi: 表中按照对应条件查找到多条记录，是否要全部更新 默认为True
        """
        try:
            self.db[collection].update(query_dict, {"$set": update_dict}, upsert=upsert, multi=multi)
        except Exception as e:
            raise Exception("更新表格{}失败！{}".format(collection, e))
        else:
            if is_show_log == "true":
                print "更新{}到表格{}成功！".format(update_dict, collection)

    def insert_one(self, collection, data_dict):
        """
        插入一条记录
        """
        collection = self.db[collection]
        return collection.insert_one(data_dict).inserted_id

    def insert_many(self, collection, data_list):
        """
        插入多条记录
        """
        collection = self.db[collection]
        collection.insert_many(data_list)

    def insert_main_table(self, collection, data):
        """
        用于插入主表
        """
        try:
            id_ = self.db[collection].insert_one(SON(data)).inserted_id
        except Exception as e:
            raise Exception("插入主表{}失败！{}".format(collection, e))
        else:
            print "插入主表{}成功！".format(collection)
        return id_
