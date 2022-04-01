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


class Mlst(Base):
    def __init__(self, bind_object):
        super(Mlst, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_mlst(self,analysis_type,desc,software,database):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "analysis_type": analysis_type,
            "status": "faile",
            "desc": desc,
            "software": software,
            "database": database,
            "params":json.dumps({"type":"mlst"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["mlst"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    @report_check
    def add_mlst_stat(self, mlst_id, file_path):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            funs = lines[0].strip().split("\t")
            for line in lines[1:]:
                lin = line.strip().split("\t")
                insert_data = {
                    "mlst_id": ObjectId(mlst_id),
                    "specimen_id": lin[0],
                    "st_id": lin[1],
                }
                for i in range(2, len(lin)):
                    insert_data[funs[i]] = lin[i]
                data_son = SON(insert_data)
                data_list.append(data_son)
        try:
            collection = self.db["mlst_stat"]
            collection.insert_many(data_list)
            main_collection = self.db["mlst"]
            main_collection.update({"_id": ObjectId(mlst_id)},
                                   {"$set": {"status": 'end','main_id': ObjectId(mlst_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_mlst_detail(self, mlst_id, file_path, specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "mlst_id": ObjectId(mlst_id),
                    "specimen_id": specimen_id,
                    "lous": line[0],
                    "ali_len": line[1],
                    "all_len": line[2],
                    "allele": line[3],
                    "gaps": line[4],
                    "coverage": line[5],
                    "identity": line[6]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["mlst_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["mlst"]
            main_collection.update({"_id": ObjectId(mlst_id)},
                                   {"$set": {"status": 'end','main_id': ObjectId(mlst_id)}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)