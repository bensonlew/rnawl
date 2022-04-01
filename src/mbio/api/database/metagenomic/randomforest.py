# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# last_modify:20180919
from biocluster.api.database.base import Base, report_check
import os,re
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
import json
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name
from bson.objectid import ObjectId

class Randomforest(Base):
    def __init__(self, bind_object):
        super(Randomforest, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_rf(self,params=None,name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'Randomforest 分析',
            'created_ts': created_ts,
            'name': name if name else "Randomforest",
            'params': params,
            'status': 'end',
        }
        try:
            collection = self.db['randomforest']
            rf_id = collection.insert_one(insert_data).inserted_id
        except Exception, e:
            #raise Exception('导入randomforest主表异常:{}'.format(e))
            self.bind_object.set_error('导入randomforest主表异常:%s', variables=(e), code="52803101")
        return rf_id

    @report_check
    def add_rf_scatter(self, rf_id,file):
        data_list = []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip('\r\n').split('\t')
                result = self.db['randomforest'].find_one({'_id': rf_id})
                task_id = result['task_id']
                samples_dic = name2id(task_id, type="task")
                insert_data = {
                    "rf_id":rf_id,
                    "specimen":samples_dic[lin[0]],
                    "x":float(lin[1]),
                    "y": float(lin[2]),
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
            try:
                collection = self.db["randomforest_scatter"]
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" % file)

    @report_check
    def add_rf_bar(self, rf_id, file):
        data_list = []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip('\r\n').split('\t')
                insert_data = {
                    "rf_id": rf_id,
                    "feature": lin[0],
                    "accuracy": float(lin[1]),
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
            try:
                collection = self.db["randomforest_bar"]
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" % file)

    @report_check
    def add_rf_evaluate(self, rf_id, file):
        data_list = []
        prob = ""
        with open(file, 'rb') as r:
            result = self.db['randomforest'].find_one({'_id': rf_id})
            task_id = result['task_id']
            samples_dic = name2id(task_id, type="task")
            head = r.readline().strip().split("\t")
            if len(head) == 2:
                data = []
                for line in r:
                    lines = line.strip().split("\t")
                    data_tmp = [lines[0], lines[1]]
                    data.append(data_tmp)
                data_list.append({"rf_id": rf_id, "table_type": "plot", "data": data})
            else:
                prob = head[2:]
                for line in r:
                    lines = line.strip().split("\t")
                    data = {
                        "rf_id": ObjectId(rf_id),
                        "table_type": "prob",
                        "name": ObjectId(samples_dic[lines[0]]),
                        "predict": lines[1],
                    }
                    for n, e in enumerate(head[2:]):
                        data[e] = lines[n + 2]
                    data_list.append(data)
        try:
            collection = self.db["randomforest_evaluate"]
            collection.insert_many(data_list)
            if prob != "":
                main_collection = self.db["randomforest"]
                main_collection.update({"_id": ObjectId(rf_id)}, {"$set": {"prob": prob,"status":'end'}})
        except Exception, e:
            self.bind_object.logger.error("导入%s信息出错:%s" % (file, e))
            self.bind_object.set_error("导入%s信息出错:%s" , variables=(file, e), code="52803102")
        else:
            self.bind_object.logger.info("导入%s信息成功!" % file)

