# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 2018.3.12

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class Crispr(Base):
    def __init__(self, bind_object):
        super(Crispr, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_crispr(self,desc):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "faile",
            "desc": desc,
            "params":json.dumps({"type":"Crispr"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["crispr"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_crispr_stat(self, crispr_id,file_path,specimen_id):
        data_list = []
        num = ''
        with open(file_path, "r") as f:
            lines = f.readlines()
            num =len(lines[1:])
        data = {
            "crispr_id": crispr_id,
            "specimen_id": specimen_id,
            "crispr_num":num
        }
        data_son = SON(data)
        data_list.append(data_son)
        try:
            collection = self.db["crispr_stat"]
            collection.insert_many(data_list)
            main_collection = self.db["crispr"]
            main_collection.update({"_id": ObjectId(crispr_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)



    @report_check
    def add_crispr_detail(self,crispr_id, file_path,specimen_id):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "crispr_id": crispr_id,
                    "specimen_id": specimen_id,
                    "location": line[0],
                    "cr_id": line[1],
                    "crispr_start": line[2],
                    "crispr_end": line[3],
                    "dr_no": int(line[4]),
                    "dr_av_len": line[5],
                    "spa_av_len": line[6]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["crispr_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["crispr"]
            main_collection.update({"_id": ObjectId(crispr_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

    @report_check
    def add_crispr_psa(self, crispr_id, file_path,specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "crispr_id": crispr_id,
                    "specimen_id": specimen_id,
                    "location": line[0],
                    "cr_id": line[1],
                    "position": line[2],
                    "dr_seq": line[3],
                    "dr_len": line[4],
                    "spa_seq": line[5],
                    "spa_len": line[6]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["crispr_psa"]
            collection.insert_many(data_list)
            main_collection = self.db["crispr"]
            main_collection.update({"_id": ObjectId(crispr_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)