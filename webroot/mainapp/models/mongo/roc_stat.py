# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'

from bson.objectid import ObjectId
import datetime
from types import StringTypes
from .core.base import Base
# from mainapp.config.db import get_mongo_client
# from biocluster.config import Config


class RocStat(Base):
    def __init__(self, bind_object):
        super(RocStat, self).__init__(bind_object)
        self._project_type = 'meta'
        # self.client = get_mongo_client()
        # self.db = self.client[Config().MONGODB]


    def create_roc(self, params, group_id=0, from_otu_table=0, level_id=0, name=None, method_type=None, top_n_id=0):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                raise Exception("from_otu_table必须为ObjectId对象或其对应的字符串!")
        if group_id != 0 and not isinstance(group_id, ObjectId):
            if isinstance(group_id, StringTypes):
                group_id = ObjectId(group_id)
            else:
                raise Exception("group_detail必须为ObjectId对象或其对应的字符串!")
        if int(level_id) not in range(1, 10):
            raise Exception("level参数%s为不在允许范围内!" % level_id)

        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        project_sn = result['project_sn']
        task_id = result['task_id']
        desc = "roc分析"
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "otu_id": from_otu_table,
            "name": name if name else "roc_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            "params": params,
            "desc": desc,
            "status": "end",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_roc"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id



