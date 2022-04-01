# -*- coding: utf-8 -*-
from __future__ import print_function
from bson.objectid import ObjectId
import types
from bson import SON
from mainapp.models.workflow import Workflow
from .core.base import Base
import random
import json
from collections import OrderedDict
import os
from biocluster.config import Config


class DeleteDemo(Base):
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(DeleteDemo, self).__init__(bind_object=self._bind_object)
        self._project_type = 'project'

    def update_main_table(self, collection_name, task_id, main_id, status):
        """
        更新主表字段
        :param collection_name:
        :param data:
        :return: 插入记录的id
        """
        conn = self.db[collection_name]
        main_id = ObjectId(main_id)
        try:
            result = conn.find_one({'_id': main_id})
            if result['task_id'] != task_id:
                info = {"success": False, "info": "The deletion condition task_id is not meet!"}
                return info
            # 同时删除demo和正式项目
            # elif result['is_get_demo'] != 'y':
            #     info = {"success": False, "info": "The deletion condition is_get_demo is not meet!"}
            #     return info
            elif result['del_mongo_status'] != "start":
                info = {"success": False, "info": "The deletion condition del_mongo_status is not meet!"}
                return info
            elif result['del_mysql_status'] != "end":
                info = {"success": False, "info": "The deletion condition del_mysql_status is not meet!"}
                return info
            else:
                db_version = result['db_version']
                conn.update({'_id': main_id}, {"$set": {'del_mongo_status': status}}, upsert=True)
                info = {"success": True, "db_version": db_version}
                return info
        except Exception as e:
            info = {"success": False, "info": "找不到数据记录!"}
            return info


if __name__ == "__main__":
    print('no test Now')
