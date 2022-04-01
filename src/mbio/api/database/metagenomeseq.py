# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import gridfs
import datetime
import os
# from biocluster.config import Config


class Metagenomeseq(Base):
    def __init__(self, bind_object):
        super(Metagenomeseq, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def add_metagenomeseq_diff(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False):
        self.bind_object.logger.info('start insert mongo zhangpeng')
        if major:
            table_id = self.create_metagenomeseq(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="51004101")
        data_list = []
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    data = [("metagenomeseq_id", table_id), ("Taxa", line_data[0]), ("samples_in_group_0", line_data[1]), ("samples_in_group_1", line_data[3]), ("counts_in_group_0", line_data[2]), ("counts_in_group_1", line_data[4]), ("oddsRatio", line_data[5]), ("lower",line_data[6]), ("upper",line_data[7]), ("fidherP",line_data[8]), ("fisherAdjP",line_data[9]), ("Pvalue",line_data[10]),("adjPvalues",line_data[11])]
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_metagenomeseq_table"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_metagenomeseq"]
            #main_collection.update({"_id":table_id}, {"$set": {"main_id": table_id}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id

    def add_metagenomeseq_list(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False):
        self.bind_object.logger.info('start insert mongo zhangpeng')
        if major:
            table_id = self.create_metagenomeseq(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="51004102")
        data_list = []
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    data = [("metagenomeseq_id", table_id), ("samples_in_group_0", line_data[0]), ("samples_in_group_1", line_data[1]), ("counts_in_group_0", line_data[2]), ("counts_in_group_1", line_data[3])]
                    data_son = SON(data)
                    data_list.append(data_son)

        try:
            conllection = self.db["sg_metagenomeseq_group_list"]
            conllection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功！"% file_path)
        return data_list, table_id








    #@report_check
    def create_metagenomeseq(self, params, group_id=0, from_otu_table=0, name=None, level_id=0):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code=51004103)
        if group_id != 0 and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                self.bind_object.set_error("group_detail必须为ObjectId对象或其对应的字符串!", code="51004104")
        if level_id not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level_id)
            self.bind_object.set_error("level参数不在允许范围内", code="51004105")

        collection = self.db["sg_otu"]  #我是不是也可以 用这个表
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "metagenomeseq分析"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else "metagenomeseq_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "params": params,
            "level_id": level_id,
            "desc": desc,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_meta_metagenomeseq"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id
