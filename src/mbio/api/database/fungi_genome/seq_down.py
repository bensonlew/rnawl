# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 2018.10.09

from biocluster.api.database.base import Base, report_check
import datetime
import json, re, os
from bson.son import SON
import pandas as pd
from collections import Counter
from bson.objectid import ObjectId


class SeqDown(Base):
    def __init__(self, bind_object):
        super(SeqDown, self).__init__(bind_object)
        self._project_type = "fungigenome"

    @report_check
    def add_seq_down(self,params=None,main_id =None,link_path=None):
        if main_id is None:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                "project_sn": project_sn,
                "task_id": task_id,
                "status": "end",
                "desc": "Seq Down",
                "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
                "created_ts": created_ts,
                "name": "序列下载",
            }
            collection = self.db["seq_down"]
            inserted_id = collection.insert_one(insert_data).inserted_id
            return inserted_id
        else:
            collection = self.db["seq_down"]
            collection.update({'_id': ObjectId(main_id)},{'$set':{'result_path':link_path,'status':'end','desc':'完成下载'}})