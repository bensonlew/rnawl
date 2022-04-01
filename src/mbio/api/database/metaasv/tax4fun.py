#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Author  : qingchen.zhang

from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
import types
import json


class Tax4fun(Base):
    def __init__(self, bind_object):
        super(Tax4fun, self).__init__(bind_object)
        self._project_type = 'metaasv'

    def add_tax4fun(self, table_name, file_path, out_type, main_id):
        main_collection_id = self.check_id(main_id)
        with open(file_path, 'rb') as r:
            lines = r.readlines()
            head = re.split('\t', lines[0].strip())
            if out_type == "level1":
                head[0] = "level1"
                specimen_list = ",".join(head[1:])
            elif out_type == "level2":
                head[0] = "level1"
                head[1] = "level2"
                specimen_list = ",".join(head[2:])
            elif out_type == "level3":
                head[0] = "level1"
                head[1] = "level2"
                head[2] = "ko"
                head[3] = "level3"
                specimen_list = ",".join(head[4:])
            else:
                head[0] = "name"
                head[1] = "description"
                specimen_list = ",".join(head[2:])
            # 更新主表specimen_list
            if out_type == "level1":
                try:
                    main_collection = self.db["tax4fun"]
                    main_collection.update({"_id": main_collection_id}, {"$set": {"specimen_list": specimen_list}})
                except Exception as e:
                    self.bind_object.logger.error("更新tax4fun表格信息{}出错:{}".format(specimen_list, e))
                else:
                    self.bind_object.logger.info("导入tax4fun表格成功:{}".format(specimen_list))
            insert_data = []
            data_temp = {
                    'type': out_type,
                    'tax4fun_id': main_collection_id
            }
            for line in lines[1:]:
                values = line.rstrip().split('\t')
                values_dict = dict(zip(head, values[0:]))
                insert_data.append(dict(data_temp, **values_dict))
        if len(insert_data) != 0 :
            try:
                collection = self.db[str(table_name)]
                collection.insert_many(insert_data)
            except Exception as e:
                self.bind_object.logger.error("导入{}表格信息出错:{}".format(str(table_name), e))
            else:
                self.bind_object.logger.info("导入{}表格成功".format(str(table_name)))
            if out_type == "level1":
                try:
                    main_collection = self.db["tax4fun"]
                    settled_params = {"software" : "R-3.3.1 (Tax4Fun)"}
                    new_specimen_list = specimen_list.split(",")
                    settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                    kegg_table_data = {"table_data": ["name", "description"] + new_specimen_list, "condition": {"type": ["KO", "enzyme"]}}
                    kegg_table_data_json = json.dumps(kegg_table_data, sort_keys=True, separators=(',', ':'))
                    level1_table_data = {"table_data": ["level1"] + new_specimen_list, "condition": {"type": "level1"}}
                    level1_table_data_json = json.dumps(level1_table_data, sort_keys=True, separators=(',', ':'))
                    level2_table_data = {"table_data": ["level1", "level2"] + new_specimen_list, "condition": {"type": "level2"}}
                    level2_table_data_json = json.dumps(level2_table_data, sort_keys=True, separators=(',', ':'))
                    level3_table_data = {"table_data": ["level1", "level2", "ko", "level3"] + new_specimen_list, "condition": {"type": "level3"}}
                    level3_table_data_json = json.dumps(level3_table_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update({"_id": main_collection_id}, {"$set": {"main_id": main_collection_id,
                                                                       "settled_params": settled_params_json,
                                                                       "kegg_table_data": kegg_table_data_json,
                                                                       "level1_table_data": level1_table_data_json,
                                                                       "level2_table_data": level2_table_data_json,
                                                                       "level3_table_data": level3_table_data_json}})
                except Exception as e:
                    self.bind_object.logger.error("更新{}表格信息出错:{}".format(str("tax4fun"), e))
        else:
            self.bind_object.logger.info("导入{}的数据为空".format(str(table_name)))

    def get_database(self, task_id):
        collection = self.db["sg_task"]
        task_info = collection.find_one({"task_id": task_id})
        if task_info:
            return task_info
        else:
            self.bind_object.logger.error("任务不存在:{}".format(str(task_id)))

    @report_check
    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！')
        return object_id
