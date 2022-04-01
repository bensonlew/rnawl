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


class Prephage(Base):
    def __init__(self, bind_object):
        super(Prephage, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_prephage(self,desc, sofware):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "faile",
            "desc": desc,
            "params":json.dumps({"type":"Prephage"}) ,
            "version": "3.1",
            "samples": "",
            "sofware": sofware,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["prephage"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        return inserted_id

    @report_check
    def add_prephage_stat(self, prephage_id, file_path, specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "prephage_id": prephage_id,
                    "specimen_id": specimen_id,
                    "len": line[1],
                    "ph_num": line[0]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["prephage_stat"]
            if len(data_list) ==1:
                collection.insert_one(data_list[0])
            else:
                collection.insert_many(data_list)
            main_collection = self.db["prephage"]
            main_collection.update({"_id": ObjectId(prephage_id)},
                                   {"$set": {"status": 'end'}})
            info = main_collection.find_one({"_id": ObjectId(prephage_id)})
            if not info:
                self.bind_object.logger.info('not found id {}'.format(prephage_id))
            else:
                self.bind_object.logger.info(' found sample {}'.format(info['samples']))
            samples = info["samples"].split(",")
            samples.append(specimen_id)
            new_samples = samples
            self.bind_object.logger.info(new_samples)
            new = ",".join(new_samples)
            self.bind_object.logger.info(new)
            main_collection.update({"_id": ObjectId(prephage_id)},{"$set": {"samples": new}},multi=True)
        except Exception, e:
            self.bind_object.set_error("prephage导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("prephage导入%s结果表成功!" % file_path)



    @report_check
    def add_prephage_detail(self,prephage_id, file_path,specimen_id):
        data_list=[]
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip('\n').split("\t")
                data = {
                    "prephage_id": prephage_id,
                    "specimen_id": specimen_id,
                    "location": line[1],
                    "ph_id": line[0],
                    "ph_start": line[2],
                    "ph_end": line[3],
                    "len": line[4],
                    "cds_num": line[7],
                    "taxon": line[5],
                    "gc": line[6],
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["prephage_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["prephage"]
            main_collection.update({"_id": ObjectId(prephage_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)

    @report_check
    def add_prephage_gene(self, prephage_id, file_path,specimen_id):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                data = {
                    "prephage_id": prephage_id,
                    "specimen_id": specimen_id,
                    "location": line[0],
                    "ph_id": line[1],
                    "gene_id": line[2],
                    "gene_start": line[4],
                    "gene_end": line[5],
                    "strand": line[3],
                    "gene_name": line[6],
                    "cog_type":line[8],
                    "nr_des": line[7]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["prephage_gene"]
            collection.insert_many(data_list)
            main_collection = self.db["prephage"]
            main_collection.update({"_id": ObjectId(prephage_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)