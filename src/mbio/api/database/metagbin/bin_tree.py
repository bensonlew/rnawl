# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.api.database.base import Base, report_check
import datetime
import json,os,re
from bson.son import SON
from bson.objectid import ObjectId

class BinTree(Base):
    def __init__(self, bind_object):
        super(BinTree, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_tree(self, params=None, name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Plot_Tree",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Plot_Tree"
        }
        collection = self.db["tree"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_tree_detail(self, insert_id, file_path, names):
        data_list = []
        with open(file_path, "r") as f:
            line = f.readlines()
            line = line[0]
            data = {
            "tree_id": ObjectId(insert_id),
            "value": line,
            }
        data_son = SON(data)
        data_list.append(data_son)
        try:
            collection = self.db["tree_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["tree"]
            main_collection.update({"_id": ObjectId(insert_id)},
                                   {"$set": {"status": 'end', "specimen": names}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)