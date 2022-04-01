# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.api.database.base import Base, report_check
import datetime
import json,re
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId

class MetagbinTree(Base):
    def __init__(self, bind_object):
        super(MetagbinTree, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_tree(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "进化树分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "metagbin tree"
        }
        collection = self.db["tree"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set': {'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_tree_detail(self, inserted_id, file, name):
        data_list = []
        with open(file, "r") as f:
            line = f.readline()
            line = re.sub('[)]([\[\]\.\;\'\"\ 0-9a-zA-Z_-]+?):', "):", line)
        data = {
            "tree_id": ObjectId(inserted_id),
            "value": line,
        }
        data_son = SON(data)
        data_list.append(data_son)
        try:
            collection = self.db["tree_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["tree"]
            main_collection.update({"_id": ObjectId(inserted_id)},
                                   {"$set": {"specimen": name}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file)