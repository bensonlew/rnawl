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


class Island(Base):
    def __init__(self, bind_object):
        super(Island, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_island(self,desc):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "failed",
            "desc": desc,
            "params":json.dumps({"type":"island"}) ,
            "version": "3.0",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "samples": ""

        }
        collection = self.db["island"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_island_stat(self, island_id, file_path, specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            gi_num = 0
            for line in lines[1:]:
                line = line.strip().split("\t")
                gi_num += int(line[1])

            data = {
                "island_id": island_id,
                "specimen_id": specimen_id,
                "gi_num": gi_num
            }
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection = self.db["island_stat"]
            collection.insert_many(data_list)
            main_collection = self.db["island"]
            info = main_collection.find_one({"_id": ObjectId(island_id)})
            samples = info["samples"].split(",")
            samples.append(specimen_id)
            new_samples = samples
            new = ",".join(new_samples)
            main_collection.update({"_id": ObjectId(island_id)},{"$set": {"samples": new}},multi=True)
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)

    @report_check
    def add_island_detail(self,island_id, file_path,specimen_id):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line_ori = line
                line = line.strip().split("\t")
                if len(line) <7:
                    self.bind_object.logger.info("缺失7列： %s" % line_ori )
                    continue
                data = {
                    "island_id": island_id,
                    "specimen_id": specimen_id,
                    "location": line[0],
                    "gi_id": line[1],
                    "island_start": line[2],
                    "island_end": line[3],
                    "len": line[4],
                    "method": line[5],
                    "cds_num": line[6],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["island_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["island"]
            main_collection.update({"_id": ObjectId(island_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

    @report_check
    def add_island_gene(self, island_id, file_path,specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "island_id": island_id,
                    "specimen_id": specimen_id,
                    "location": line[0],
                    "gi_id": line[1],
                    "gene_id": line[2],
                    "gene_start": line[3],
                    "gene_end": line[4],
                    "strand": line[5],
                    "product": line[6],
                    "name": line[7],
                    "cog_type": line[8]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["island_gene"]
            collection.insert_many(data_list)
            main_collection = self.db["island"]
            main_collection.update({"_id": ObjectId(island_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)