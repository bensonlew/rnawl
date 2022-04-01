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


class Plasmid(Base):
    def __init__(self, bind_object):
        super(Plasmid, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_plasmid(self,desc,software,database):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "faile",
            "desc": desc,
            "software": software,
            "database": database,
            "params":json.dumps({"type":"Plasmid"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["plasmid"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    @report_check
    def add_plasmid_stat(self, plasmid_id, file_path):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                num = int(line[1])+int(line[2])+int(line[3])
                data = {
                    "plasmid_id": plasmid_id,
                    "specimen_id": line[0],
                    "seq_num":num,
                    "chr_num": line[2],
                    "known_pla_num": line[1],
                    "new_pla_num": line[3],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["plasmid_stat"]
            collection.insert_many(data_list)
            main_collection = self.db["plasmid"]
            main_collection.update({"_id": ObjectId(plasmid_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_plasmid_detail(self, plasmid_id, file_path, dict):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                location = dict[line[1]]
                data = {
                    "specimen_id": line[0],
                    "plasmid_id": plasmid_id,
                    "location": location,
                    "type": line[2],
                    "phylum": line[3],
                    "probability": line[4],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["plasmid_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["plasmid"]
            main_collection.update({"_id": ObjectId(plasmid_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_plasmid_anno(self, plasmid_id, file_path,specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "plasmid_id": plasmid_id,
                    "specimen_id": specimen_id,
                    "location": line[0],
                    "seq_id": line[1],
                    "gc": line[2],
                    "seq_len": line[3],
                    "identity": line[4],
                    "evalue": line[5],
                    "acc_num": line[6],
                    "plasmid_type": line[8],
                    "taxon": line[7],
                    "strain": line[9],
                    "plsa_name": line[10],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["plasmid_anno"]
            collection.insert_many(data_list)
            main_collection = self.db["plasmid"]
            main_collection.update({"_id": ObjectId(plasmid_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)