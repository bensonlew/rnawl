# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId

class TetraPcoa(Base):
    def __init__(self, bind_object):
        super(TetraPcoa, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_pcoa(self, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "tetra pcoa 分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": "bin_tetra"
        }
        collection = self.db["tnf"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set': {'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_pcoa_detail(self, inserted_id, file, file2):
        data_list = []
        with open (file,'r') as f:
            lines = f.readlines()
            for lin in lines[1:]:
                line =lin.rstrip('\r\n').split('\t')
                data ={
                    "pcoa_id": ObjectId(inserted_id),
                    "scaf_id": line[0],
                    "pc1": line[1],
                    "pc2": line[2],
                    "pc3": line[3],
                }
                data_son = SON(data)
                data_list.append(data_son)
        with open (file2, 'r') as f:
            lines = f.readlines()
            pc1 =lines[1].rstrip('\r\n').split('\t')[1]
            pc2 = lines[2].rstrip('\r\n').split('\t')[1]
            pc3 = lines[3].rstrip('\r\n').split('\t')[1]
            importance = {"pc1": pc1, "pc2": pc2, "pc3": pc3}
        try:
            collection = self.db["tnf_detail"]
            collection.insert_many(data_list)
            collection_db = self.db["tnf"]
            collection_db.update({'_id': ObjectId(inserted_id)}, {'$set': {'importance': importance}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file)