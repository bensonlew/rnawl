# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId

class InquireSeq(Base):
    def __init__(self, bind_object):
        super(InquireSeq, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_seq(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Inquire_seq",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "Inquire_Seq"
        }
        collection = self.db["seq"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': ObjectId(inserted_id)},{'$set':{'main_id': ObjectId(inserted_id)}})
        return inserted_id

    @report_check
    def add_anno_nr_detail(self, inserted_id, anno):
        data_list = []
        with open (anno,'r') as f:
            lines = f.readlines()
            for lin in lines:
                line =lin.rstrip('\r\n').split('\t')
                data ={
                    "seq_id": ObjectId(inserted_id),
                    "scf_id": line[0],
                    "start": line[6],
                    "end": line[7],
                    "identity": line[2],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["seq_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (anno, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % anno)