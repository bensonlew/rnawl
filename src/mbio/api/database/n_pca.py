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


class NPca(Base):
    def __init__(self, bind_object):
        super(NPca, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def add_environmental_regression_site(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False):
        self.bind_object.logger.info('start insert mongo zhangpeng')
        if major:
            table_id = self.create_environmental_regression(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！", code="51004301")
        data_list = []
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    data = [("sample_name", line_data[0]), ("X_PCA", line_data[1]), ("Y_factor", line_data[2])]
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_environmental_regression_curve"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id

    # @report_check
    def add_n_pca_site(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False):
        if major:
            table_id = self.create_environmental_regression(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51004302")
        data_list = []
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[0:]
            data_temp = []
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                #insert_data = {
                    #'n_pca_id':updata_id,
                    #'type':table_type
                #}
                #insert_data['name']=values[0]
                #values_dict = dict(zip(columns,values[1:]))
                #a<-len(values)

                data_temp = zip(columns,values[1:])
                #print data_temp
                #collection.insert_many(data_temp)
                data_son = SON(data_temp)
                data_list.append(data_son)
        try:
            collection = self.db["sg_n_pca_sites"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错", code="51004303")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list




    # @report_check
    def add_n_pca_max(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False):
        if major:
            table_id = self.create_environmental_regression(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51004304")
        data_list = []
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[0:]
            data_temp = []
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                #insert_data = {
                    #'n_pca_id':updata_id,
                    #'type':table_type
                #}
                #insert_data['name']=values[0]
                #values_dict = dict(zip(columns,values[1:]))
                #a<-len(values)

                data_temp = zip(columns,values[1:])
                #print data_temp
                #collection.insert_many(data_temp)
                data_son = SON(data_temp)
                data_list.append(data_son)
        try:
            collection = self.db["sg_n_pca_sites_max"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list


    # @report_check
    def add_n_pca_min(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False):
        if major:
            table_id = self.create_environmental_regression(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51004305")
        data_list = []
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[0:]
            data_temp = []
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                #insert_data = {
                    #'n_pca_id':updata_id,
                    #'type':table_type
                #}
                #insert_data['name']=values[0]
                #values_dict = dict(zip(columns,values[1:]))
                #a<-len(values)

                data_temp = zip(columns,values[1:])
                #print data_temp
                #collection.insert_many(data_temp)
                data_son = SON(data_temp)
                data_list.append(data_son)
        try:
            collection = self.db["sg_n_pca_sites_min"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错", code="51004306")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list


    # @report_check
    def add_n_pca_importance(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False):
        if major:
            table_id = self.create_environmental_regression(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51004307")
        data_list = []
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[0:]
            #data_temp = []
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                #insert_data = {
                    #'n_pca_id':updata_id,
                    #'type':table_type
                #}
                #insert_data['name']=values[0]
                #values_dict = dict(zip(columns,values[1:]))
                #a<-len(values)
                columns1 = []
                for ka in columns:
                    columns1.append(ka.strip("\""))
                main_collection = self.db["sg_npca"]
                main_collection.update({"_id": table_id}, {"$set":{"columns": columns1}})
                data_temp = zip(columns1,values[1:])
                a_id = ("npca_id", table_id)
                data_temp.append(a_id)
                #print data_temp
                #collection.insert_many(data_temp)
                data_son = SON(data_temp)
                data_list.append(data_son)
        try:
            collection = self.db["sg_npca_importance"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错", code="51004306")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    # @report_check
    def add_n_pca_sitesall(self, file_path, table_id = None, group_id = None, from_otu_table = None, level_id = None, major = False):
        if major:
            table_id = self.create_environmental_regression(self, params, group_id, from_otu_table, level_id)
        else:
            #if table_id is None:
                #raise Exception("major为False时需提供table_id!")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!", code="51004307")
        data_list = []

        a1 = len(file_path)
        a2 = len(file_path[0])
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    data = [("npca_id", table_id), ("newnames", line_data[1].strip("\"")), ("sample", line_data[2].strip("\"")), ("highgroup", line_data[3].strip("\"")), ("lowgroup", line_data[4].strip("\"")), ("pca_value", line_data[5]), ("pcai", line_data[6].strip("\"")), ("pca_max_value", line_data[7]), ("pca_min_value", line_data[8])]
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["sg_npca_site"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_npca"]
            #main_collection.update({"_id":table_id}, {"$set":{"main_id": table_id}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错", code="51004306")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list


    #@report_check
    def create_environmental_regression(self, params, group_id=0, from_otu_table=0, name=None, level_id=0):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51004308")
        if group_id != 0 and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                self.bind_object.set_error("group_detail必须为ObjectId对象或其对应的字符串!", code="51004309")
        if level_id not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level_id)
            self.bind_object.set_error("level参数不在范围内", code="51004310")

        collection = self.db["sg_otu"]  #我是不是也可以 用这个表
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "roc分析"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": name if name else "oturoc_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "params": params,
            "level_id": level_id,
            "desc": desc,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_meta_roc"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id
