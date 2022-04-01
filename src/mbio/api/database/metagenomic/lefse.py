# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last modify: 2017.11.17

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
from mainapp.libs.param_pack import group_detail_sort


class Lefse(Base):
    def __init__(self, bind_object):
        super(Lefse, self).__init__(bind_object)
        self._project_type = "metagenomic"


    @report_check
    def add_lefse(self, anno_type, geneset_id, anno_id, level_id,  group_id, name=None, params=None):
        task_id = self.bind_object.sheet.id
        group_detail = {'All': [str(i) for i in spname_spid.values()]}
        params['group_detail'] = group_detail_sort(group_detail)
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "anno_type": anno_type,
            "geneset_id": ObjectId(geneset_id),
            'group_id': group_id,
            "anno_id": anno_id,
            "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "level_id": level_id,
            "status": "faile",
            "desc": "lefse分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["lefse"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_lefse_detail(self, file_path, lefse_id):
        data_list = []
        with open(file_path, "r") as f:
            envs = f.readline().strip().split("\t")[1:]
            for line in f:
                line = line.strip().split("\t")
                if line[3] == '':
                    lda = 0
                else:
                    lda = float(line[3])
                if line[4] != '-':
                    line[4] = float(line[4])
                data = {
                    "species_lefse_id": ObjectId(lefse_id),
                    "species_name": line[0],
                    "category_name": line[2],
                    "median": float(line[1]),
                    "lda": lda,
                    "pvalue": line[4]
                }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["lefse_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["lefse"]
            main_collection.update({"_id": ObjectId(lefse_id)},
                                   {"$set": {"status": 'end'}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
            self.bind_object.set_error("导入lefse结果表出错", code="52801401")
        else:
            self.bind_object.logger.info("导入%s结果表成功!" %file_path)
