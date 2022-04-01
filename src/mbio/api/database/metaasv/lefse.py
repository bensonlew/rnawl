# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
import datetime
import json


class Lefse(Base):
    def __init__(self, bind_object):
        super(Lefse, self).__init__(bind_object)
        self._project_type = "metaasv"

    @report_check
    def add_species_difference_lefse_detail(self, file_path, params=None, group_id=None, from_otu_table=None, table_id=None, major=False):
        if major:
            table_id = self.create_species_difference_lefse(params, group_id, from_otu_table)
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
            r.readline()
            is_use = False
            for line in r:
                line = line.strip('\n')
                line_data = line.split('\t')
                data = [
                    ("lefse_id", table_id),
                    ("species_name", line_data[0]),
                    ("specimen", line_data[2]),
                    ("median", float(line_data[1])),
                ]
                if line_data[3] != "":
                    data.append(("lda", float(line_data[3])))
                if line_data[4] not in ["-", ""]:
                    data.append(("pvalue", float(line_data[4])))
                data_son = SON(data)
                data_list.append(data_son)
                if not line_data[2] == '-':
                    is_use = True
        if is_use:
            coll_main = self.db["lefse"]
            coll_main.update({"_id": ObjectId(table_id)}, {"$set": {"lda_png_id": 'is_useful',
                                                                    "main_id": table_id,
                                                                    "lda_cladogram_id": 'is_useful'}})
        try:
            collection = self.db["lefse_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        try:
            main_collection = self.db["lefse"]
            settled_params = {"software" : "R-3.3.1"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            table_data = {
                "table_data": ["species_name","specimen","median", "lda", "pvalue"]}
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
            column_data = {
                "column_data": {"name":"species_name",
                                "data":"median",
                            "category": ""}}
            column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": table_id},{"$set": {"main_id": table_id,
                                                              "settled_params": settled_params_json,
                                                              "table_data": table_data_json,
                                                              "column_data": column_data_json}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入信息出错")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file_path)
        return data_list, table_id


    @report_check
    def create_species_difference_lefse(self, params, group_id=0, from_otu_table=0, name=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        if group_id != 0 and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "lefse分析"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "name": name if name else "LEfSe_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "params": params,
            "lda_cladogram_id": "",
            "lda_png_id": "",
            "desc": desc,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["lefse"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id
