# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20181214


from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
from types import StringTypes
from biocluster.config import Config
from bson.son import SON

class AssessGc(Base):
    def __init__(self, bind_object):
        super(AssessGc, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_assembly_gc(self,genome_id,gc_path=None, params=None, name=None):
        old_task_id = str(self.bind_object.sheet.id).split('_')[0]
        num = str(self.bind_object.sheet.id).split('_')[1]
        task_id = old_task_id + "_" + num
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Gc_depth统计主表",
            "params":json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "genome_id":genome_id,# G_bin_01
            "name": name if name else "Gc_depth_Origin",
            "path": gc_path,
        }
        collection = self.db["assembly_gc"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id': inserted_id}})
        return inserted_id