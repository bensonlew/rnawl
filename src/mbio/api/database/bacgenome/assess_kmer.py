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


class AssessKmer(Base):
    def __init__(self, bind_object):
        super(AssessKmer, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_assess_kmer(self,specimen_id,desc):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "faile",
            "desc": desc,
            "params":json.dumps({"type":"kmer frequency"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "specimen_id": specimen_id
        }
        collection = self.db["assess_kmer"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_assess_kmer_detail(self,assess_kmer_id, file_path):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "assess_kmer_id": assess_kmer_id,
                    "depth": line[0],
                    "frequency": float(line[1])
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assess_kmer_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["assess_kmer"]
            main_collection.update({"_id": ObjectId(assess_kmer_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

