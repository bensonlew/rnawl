# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'

from biocluster.api.database.base import Base, report_check
import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes


class MetabsetVenn(Base):
    def __init__(self, bind_object):
        super(MetabsetVenn, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_venn(self, name=None, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "正在计算",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),

        }
        collection = self.db["metabset_venn"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': inserted_id}, {'$set': {'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_venn_detail(self, venn_path, venn_id):
        data_list = []
        with open(venn_path, 'rb') as r:
            line0 = r.readline()
            g_num = len(line0.split('\t')[0].split('&'))

        with open(venn_path, 'rb') as r:
                for line in r:
                    line = line.strip().split("\t")
                    c_num = len(line[0].split('&'))
                    if g_num > 6:
                        if c_num not in [1, g_num]:
                            continue
                    insert_data = {
                        'venn_id': ObjectId(venn_id),
                        'label': line[0],
                        'num': int(line[1]),
                        }
                    if int(line[1]):
                       insert_data['content_list'] = line[2]
                    else:
                       insert_data['content_list'] = ""
                    data_list.append(insert_data)
        try:
            collection = self.db["metabset_venn_detail"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入venn数据出错:%s" % e)
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
            collection = self.db["metabset_venn_graph"]
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入Venn画图数据出错:%s" % e)
        else:
            self.bind_object.logger.error("导入Venn画图数据成功")
