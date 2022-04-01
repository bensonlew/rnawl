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


class HouseKeep(Base):
    def __init__(self, bind_object):
        super(HouseKeep, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_house(self,desc,software,database):
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
        collection = self.db["house_keep"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    @report_check
    def add_house_detail(self, house_id, file_path,specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "house_id": house_id,
                    "specimen_id": specimen_id,
                    "gene_id":line[0],
                    "house_gene": line[1],
                    "location": line[2],
                    "strand": line[3],
                    "start": line[4],
                    "end": line[5],
                    "identity": line[6],
                    "coverage": line[7],
                    "evalue": line[8],
                    "score": line[9],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["house_keep_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["house_keep"]
            main_collection.update({"_id": ObjectId(house_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)
