# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

from biocluster.api.database.base import Base, report_check
import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes


class Venn(Base):
    def __init__(self, bind_object):
        super(Venn, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_venn(self, anno_type, name=None, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "anno_type": anno_type,
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
        collection = self.db["venn"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_venn_detail(self, venn_path, venn_id):
        data_list = []
        with open(venn_path, 'rb') as r:
                for line in r:
                    line = line.strip().split("\t")
                    insert_data = {
                        'venn_id': ObjectId(venn_id),
                        'label': line[0],
                        #'content_list': line[2],
                        'num': int(line[1]),
                        }
                    if int(line[1]):
                       insert_data['content_list'] = line[2]
                    else:
                       insert_data['content_list'] = ""
                    data_list.append(insert_data)
        try:
            collection = self.db["venn_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入venn数据出错:%s" % e)
            self.bind_object.set_error("导入venn数据出错", code="52803001")
        else:
            self.bind_object.logger.info("导入venn数据成功")

    @report_check
    def add_venn_graph(self, venn_graph_path, venn_id):
        data_list = []
        with open(venn_graph_path, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                insert_data = {
                    'venn_id': ObjectId(venn_id),
                    'category_name': line[0],
                    'content_list':  line[1]
                    }
                data_list.append(insert_data)
        try:
            collection = self.db["venn_graph"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入Venn画图数据出错:%s" % e)
            self.bind_object.set_error("导入Venn画图数据出错", code="52803002")
        else:
            self.bind_object.logger.error("导入Venn画图数据成功")

    @report_check
    def add_venn_pie(self, venn_pie_path, venn_id):  # venn图中画饼图的数据
        data_list = []
        with open(venn_pie_path, "r") as k:
            for line in k:
                line = line.strip().split("\t")
                #pei_data = line[2].replace("}, {",", ").lstrip('[').rstrip(']')
                insert_data = {
                    'venn_id': ObjectId(venn_id),
                    'category_name': line[0],
                    'others': line[1],
                    'venn_pie_data': eval(line[2])
                    #'venn_pie_data': eval(pei_data)
                    }
                data_list.append(insert_data)
        try:
            collection = self.db["venn_pie_data"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入Venn饼图数据出错:%s" % e)
            self.bind_object.set_error("导入Venn饼图数据出错", code="52803003")
        else:
            self.bind_object.logger.error("导入Venn饼图数据成功")
