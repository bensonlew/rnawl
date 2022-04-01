# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 20210527

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class Srna(Base):
    def __init__(self, bind_object):
        super(Srna, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_srna(self, desc,software,database):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "faile",
            "desc": desc,
            "software": software,
            "database": database,
            "params":json.dumps({"type":"srna"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["srna"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    @report_check
    def add_srna_stat(self, srna_id, file_path,specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "srna_id": ObjectId(srna_id),
                    "specimen_id": specimen_id,
                    "srna_num":line[0],
                    "total_len": line[1],
                    "percent": line[2],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["srna_stat"]
            collection.insert_many(data_list)
            main_collection = self.db["srna"]
            main_collection.update({"_id": ObjectId(srna_id)},
                                   {"$set": {"status": 'end',"main_id": ObjectId(srna_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_srna_detail(self, srna_id, file_path):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "srna_id": ObjectId(srna_id),
                    "specimen_id": line[0],
                    "location": line[1],
                    "srna": line[2],
                    "rfam_id": line[3],
                    "family": line[4],
                    "des": line[5],
                    "start": int(line[6]),
                    "end": int(line[7]),
                    "srna_len": int(line[8]),
                    "strand": line[9],
                    "evalue": float(line[10]),
                    "score": float(line[11]),
                    "seq": line[12],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["srna_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["srna"]
            main_collection.update({"_id": ObjectId(srna_id)},
                                   {"$set": {"status": 'end',"main_id":ObjectId(srna_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)