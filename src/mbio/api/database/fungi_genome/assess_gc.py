# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 2018.3.5

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class AssessGc(Base):
    def __init__(self, bind_object):
        super(AssessGc, self).__init__(bind_object)
        self._project_type = "fungigenome"


    @report_check
    def add_assess_gc(self,desc):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "faile",
            "desc": desc,
            "params":json.dumps({"type":"gc_depth"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["assess_gc"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_assess_gc_detail(self,assess_gc_id, sample_name,window,file_path):
        data_list=[]
        data = {
            "assess_gc_id": assess_gc_id,
            "specimen_id": sample_name,
            "window": window,
            "path": file_path
        }
        data_son = SON(data)
        data_list.append(data_son)
        try:
            collection = self.db["assess_gc_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["assess_gc"]
            main_collection.update({"_id": ObjectId(assess_gc_id)},
                                   {"$set": {"status": 'end',"origin_id": ObjectId(assess_gc_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)
