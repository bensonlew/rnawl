# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'


from biocluster.api.database.base import Base, report_check
import os
import json
import datetime
from bson.objectid import ObjectId


class Micropita(Base):
    def __init__(self, bind_object):
        super(Micropita, self).__init__(bind_object)
        self._project_type = "meta"

    @report_check
    def add_micropita(self, otu_id, level_id, name=None, params=None):
        task_id = self.bind_object.sheet.id
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "otu_id": ObjectId(otu_id),
            "name": name,
            "level_id": int(level_id),
            "status": "end",
            "desc": "MicroPITA分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_micropita"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_micropita_detail(self, micropita_id, pcoa_file, select_file):
        select_method = []
        data_list = []
        with open(pcoa_file, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "micropita_id": ObjectId(micropita_id),
                    "type": "plot",
                    "PC1": line[0],
                    "PC2": line[1],
                    "name": line[2],
                }
                data_list.append(data)
        with open(select_file, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "micropita_id": ObjectId(micropita_id),
                    "type": "pick",
                    "name": line[0],
                    "specimen_list": ';'.join(line[1:]),
                }
                select_method.append(line[0])
                data_list.append(data)
        try:
            collection = self.db["sg_micropita_detail"]
            collection.insert_many(data_list)
            self.bind_object.logger.info("开始刷新主表写树")
            main_collection = self.db["sg_micropita"]
            main_collection.update({"_id": ObjectId(micropita_id)}, {"$set": {"select_method": select_method}})

            os.remove(pcoa_file)
        except Exception, e:
            self.bind_object.logger.error("导入micropita详情表出错:%s" % e)
            self.bind_object.set_error("导入micropita出错", code="51004201")
        else:
            self.bind_object.logger.info("导入micropita详情表表格成功")
