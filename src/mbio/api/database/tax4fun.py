#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Author  : houshuang 2019/9/27

from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
import types


class Tax4fun(Base):
    def __init__(self, bind_object):
        super(Tax4fun, self).__init__(bind_object)
        self._project_type = 'meta'

    def add_tax4fun(self, table_name, file_path, out_type, main_id):
        main_id = self.check_id(main_id)
        with open(file_path, 'rb') as r:
            lines = r.readlines()
            head = re.split('\t', lines[0].strip())
            if out_type == "level1":
                head[0] = "level1"
                specimen_list = ",".join(head[1:])
            elif out_type == "level2":
                head[0] = "level1"
                head[1] = "level2"
            elif out_type == "level3":
                head[0] = "level1"
                head[1] = "level2"
                head[2] = "ko"
                head[3] = "level3"
            else:
                head[0] = "name"
                head[1] = "description"
            # 更新主表specimen_list
            if out_type == "level1":
                try:
                    main_collection = self.db["sg_tax4fun"]
                    main_collection.update({"_id": main_id}, {"$set": {"specimen_list": specimen_list}})
                except Exception as e:
                    self.bind_object.logger.error("更新sg_tax4fun表格信息{}出错:{}".format(specimen_list, e))
                else:
                    self.bind_object.logger.info("导入sg_tax4fun表格成功:{}".format(specimen_list))

            insert_data = []
            data_temp = {
                    'type': out_type,
                    'tax4fun_id': main_id
            }
            for line in lines[1:]:
                values = line.rstrip().split('\t')
                values_dict = dict(zip(head, values[0:]))
                insert_data.append(dict(data_temp, **values_dict))
            if insert_data:
                try:
                    collection = self.db[str(table_name)]
                    collection.insert_many(insert_data)
                    main_collection = self.db['sg_tax4fun']
                    #main_collection.update({"_id": main_id}, {"$set": {"main_id": main_id}})
                except Exception as e:
                    self.bind_object.logger.error("导入{}表格信息出错:{}".format(str(table_name), e))
                else:
                    self.bind_object.logger.info("导入{}表格成功".format(str(table_name)))
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
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', code="51008601")
        return object_id
