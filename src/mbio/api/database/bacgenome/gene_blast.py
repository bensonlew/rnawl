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


class GeneBlast(Base):
    def __init__(self, bind_object):
        super(GeneBlast, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_blast(self,params=None,name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "基因搜索",
            "params":json.dumps({"type":"GeneBlast"}) ,
            "name": name if name else "GeneBlast",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["blast"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    @report_check
    def add_blast_detail(self, blast_id, file_path, nr_id):
        data_list=[]
        nr_collection = self.db['anno_nr_detail']
        if not isinstance(blast_id, ObjectId):
            blast_id = ObjectId(blast_id)
        if not isinstance(nr_id, ObjectId):
            nr_id = ObjectId(nr_id)
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "blast_id": blast_id,
                    "specimen_id": line[0],
                    "gene_id": line[1],
                    "location": line[2],
                    "q_start": int(float(line[3])),
                    "q_end": int(float(line[4])),
                    "s_start": int(float(line[5])),
                    "s_end": int(float(line[6])),
                    "alien_length": int(float(line[7])),
                    "identity": float(line[8]),
                    "coverage": float(line[9]),
                    "score": float(line[11]),
                    "evalue": line[10],
                }
                nr_detail =  nr_collection.find_one({"nr_id":nr_id,"specimen_id": line[0], "gene_id": line[1]})
                nr_des = nr_detail["gene_des"]
                data["nr_des"] = nr_des
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["blast_detail"]
            if len(data_list) > 0:
                collection.insert_many(data_list)
            main_collection = self.db['blast']
            main_collection.update({'_id':blast_id},{"$set":{"main_id":blast_id}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

