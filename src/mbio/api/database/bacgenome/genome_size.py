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


class GenomeSize(Base):
    def __init__(self, bind_object):
        super(GenomeSize, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_assess_size(self,desc):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "faile",
            "desc": desc,
            "params":json.dumps({"type":"genome_size"}) ,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["assess_size"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_assess_size_detail(self,assess_size_id, file_path,specimen_id):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "assess_size_id": assess_size_id,
                    "specimen_id": specimen_id,
                    "method": line[0],
                    "kmer": line[1],
                    "k_indiv": line[2],
                    "k_depth": line[3],
                    "size": line[4],
                    "r_size": line[5],
                    "h_rate": line[6],
                    "repeat": line[7],
                    "same_kmer": line[8],
                    "accuracy_rate": line[9],
                    "u_base": line[10],
                    "aver_depth": line[11]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assess_size_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["assess_size"]
            main_collection.update({"_id": ObjectId(assess_size_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

