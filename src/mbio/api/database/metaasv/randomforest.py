# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import json
import datetime


class Randomforest(Base):
    def __init__(self, bind_object):
        super(Randomforest, self).__init__(bind_object)
        self._project_type = 'metaasv'

    @report_check
    def add_randomforest_error(self, file_path, table_id=None, group_id=None, from_otu_table=None, level_id=None,
                               major=False, params=None):
        self.bind_object.logger.info('start insert mongo')
        if major:
            table_id = self.create_randomforest(self, params, group_id, from_otu_table, level_id)
        else:
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        data_list1 = []
        b1 = len(file_path)
        b2 = len(file_path[0])
        with open(file_path, "rb") as r:
            for line in r:
                line = line.strip('\n')
                line_data = line.split('\t')
                data_list1.append(line_data)
            data = [("randomforest_id", table_id), ("class.error", line_data[b2 + 1])]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["randomforest"]
            collection.insert_many(data_list)

            # main_collection = self.db["randomforest"]
            #main_collection.update_one({'_id': table_id},{'$set': {'main_id': table_id}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id

    @report_check
    def add_randomforest_dim(self, file_path, table_id=None, group_id=None, from_otu_table=None, level_id=None,
                             major=False,params=None):
        if major:
            table_id = self.create_randomforest(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file_path, 'rb') as r:
            i = 0
            for line in r:
                if i == 0:
                    i = 1
                else:
                    line = line.strip('\n')
                    line_data = line.split('\t')
                    data = [("randomforest_id", table_id),
                            ("specimen_name", line_data[0]),
                            ("x", line_data[1]),
                            ("y", line_data[2])]
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["randomforest_scatter"]
            collection.insert_many(data_list)
            main_collection = self.db["randomforest"]
            settled_params = {"software" : "R-3.3.1 (randomForest)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            scatter_data = {
                "scatter_data": {"name": "specimen_name",
                "category": ""}
            }
            scatter_data_json = json.dumps(scatter_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": ObjectId(table_id)}, {"$set": {"main_id": table_id,
                                                                          "settled_params": settled_params_json,
                                                                          "scatter_data": scatter_data_json}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    @report_check
    def add_randomforest_vip(self, file_path, table_id=None, group_id=None, from_otu_table=None, level_id=None,
                             major=False,params=None):
        if major:
            table_id = self.create_randomforest(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
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
                    data = [("randomforest_id", table_id),
                            ("species_name", line_data[0]),
                            ("accuracy", float(line_data[1]))]
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["randomforest_bar"]
            collection.insert_many(data_list)
            main_collection = self.db["randomforest"]
            column_data = {
                "column_data": {"name": "specimen_name",
                "data": "accuracy"}
            }
            column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
            bar_table_data = {
                "table_data": ["species_name","accuracy" ]
            }
            bar_table_data_json = json.dumps(bar_table_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": ObjectId(table_id)}, {"$set": {"main_id": table_id,
                                                                          "bar_table_data": bar_table_data_json,
                                                                          "column_data": column_data_json}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list

    @report_check
    def add_randomforest_evaluate(self, file_path, table_id=None, group_id=None, from_otu_table=None, level_id=None,
                                  major=False,params=None):
        if major:
            table_id = self.create_randomforest(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        with open(file_path, 'rb') as r:
            head = r.readline().strip().split("\t")
            if len(head) == 2:## 选择AUC值等方法得到的结果
                for line in r:
                    line_split = line.strip().split("\t")
                    data_list.append({"randomforest_id": table_id,
                                      "name":"",
                                      "x": int(line_split[0]),
                                      "y": float(line_split[1])})
        try:
            collection = self.db["randomforest_line_chart"]
            collection.insert_many(data_list)
            main_collection = self.db["randomforest"]
            line_data = {
                "line_data": {"name": "name",
                "category": ""}
            }
            line_data_json = json.dumps(line_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": ObjectId(table_id)}, {"$set": {"main_id": table_id,
                                                                          "line_data": line_data_json}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)

    def add_randomforest_predict(self, file_path, table_id=None, group_id=None, from_otu_table=None, level_id=None,
                                  major=False,params=None):
        if major:
            table_id = self.create_randomforest(self, params, group_id, from_otu_table, level_id)
        else:
            if table_id is None:
                self.bind_object.set_error("major为False时需提供table_id!")
            if not isinstance(table_id, ObjectId):
                if isinstance(table_id, StringTypes):
                    table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("table_id必须为ObjectId对象或其对应的字符串!")
        data_list = []
        prob = ""
        with open(file_path, 'rb') as r:
            head = r.readline().strip().split("\t")
            prob = head[2:]
            for line in r:
                lines = line.strip().split("\t")
                data = {
                    "randomforest_id": ObjectId(table_id),
                    "name": lines[0],
                    "predict": lines[1],
                }
                for n, e in enumerate(head[2:]):
                    data[e] = lines[n + 2]
                data_list.append(data)
        try:
            collection = self.db["randomforest_table"]
            collection.insert_many(data_list)
            if prob != "":
                main_collection = self.db["randomforest"]
                table_data = {"table_data": ["name", "predict"] + prob}
                table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": ObjectId(table_id)}, {"$set": {"prob": prob,
                                                                              "main_id": table_id,
                                                                              "table_data": table_data_json}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)


    @report_check
    def create_randomforest(self, params, group_id=0, from_otu_table=0, name=None, level_id=0):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        if group_id != 0 and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                self.bind_object.set_error("group_detail必须为ObjectId对象或其对应的字符串!")
        if level_id not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level_id)
            self.bind_object.set_error("level不在允许范围内")

        collection = self.db["asv"]  # 我是不是也可以 用这个表
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "randomforest分析"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": name if name else "randomforest_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "params": params,
            "level_id": level_id,
            "desc": desc,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["randomforest"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id
