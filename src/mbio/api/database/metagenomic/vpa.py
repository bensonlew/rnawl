# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'

from biocluster.api.database.base import Base, report_check
import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes


class Vpa(Base):
    def __init__(self, bind_object):
        super(Vpa, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_vpa(self, name=None, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            #"anno_type": anno_type,
            #"geneset_id": geneset_id,
            #"anno_id": anno_id,
            #"level_id": level_id,
            #"group_id": group_id,
            "status": "end",
            "desc": "正在计算",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["vpa"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_vpa_detail(self, vpa_path, vpa_id):
        data_list = []
        if not isinstance(vpa_id, ObjectId):
            vpa_id = ObjectId(vpa_id)
        with open(vpa_path, 'rb') as r:
                r.readline()
                for line in r:
                    line = line.strip().split("\t")
                    insert_data = {
                        'vpa_id': vpa_id,
                        'env' : line[0],
                        'radj': line[1],
                        'group' : line[2]
                        }
                    data_list.append(insert_data)
        try:
            collection = self.db["vpa_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入Vpa数据出错:%s" % e)
            self.bind_object.set_error("import Vpa data ERROR", code="52804201")
        else:
            self.bind_object.logger.info("导入Vpa数据成功")

    @report_check
    def add_vpa_graph(self, vpa_graph_path, vpa_id, png_path):
        data_list = []
        if not isinstance(vpa_id, ObjectId):
            vpa_id = ObjectId(vpa_id)
        with open(vpa_graph_path, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                insert_data = {
                    'vpa_id': vpa_id,
                    'label': line[0].strip('"'),
                    'num' : line[1]
                    }
                data_list.append(insert_data)
        try:
            collection = self.db["vpa_graph"]
            collection.insert_many(data_list[:-1])
            main=self.db['vpa']
            main.update_one({'_id':vpa_id},{"$set":{'residual':data_list[-1]['num'], 'pic':png_path}})

        except Exception, e:
            self.bind_object.logger.error("导入Vpa画图数据出错:%s" % e)
            self.bind_object.set_error("import Vpa graph data ERROR", code="52804202")
        else:
            self.bind_object.logger.info("导入Vpa画图数据成功")

