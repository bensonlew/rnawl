# -*- coding: utf-8 -*-
# __author__ = 'liulinmeng'
# last_modify 20181115

import os,re
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check
#from biocluster.config import Config

class Roc(Base):
    def __init__(self, bind_object):
        super(Roc, self).__init__(bind_object)
        self._project_type = "meta"

    @report_check
    def add_roc(self, task_id, project_sn, params = None, name = None):
        task_id = task_id
        project_sn = project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'Roc analysis',
            'created_ts': created_ts,
            'name': name if name else "ROC",
            'params': params,
            'status': 'start'
        }
        try:
            collection = self.db['sg_roc']
            roc_id = collection.insert_one(insert_data).inserted_id
        except Exception, e:
            #raise Exception('Error when import ROC main_table:{}'.format(e))
            self.bind_object.logger.info("导入ROC主表出错: %s" % e)
        else:
            self.bind_object.logger.info("导入ROC主表成功!")
            return roc_id

    @report_check
    def add_roc_curve(self, roc_id, type, file):
        data_list = []
        with open(file, "r") as f:
            lines = f.readlines()
            n = 1
            for line in lines[1:]:
                lst = line.strip("\r\n").split("\t")
                insert_data = {
                    "x": (100-float(lst[0]))/100.0,
                    "y": float(lst[1])/100.0,
                    "type": type,
                    "roc_id": ObjectId(roc_id),
                    "sort": n,
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
                n += 1
            try:
                collection = self.db['sg_roc_curve']
                collection.insert_many(data_list)

                main_collection = self.db["sg_roc"]
                #main_collection.update({"_id": ObjectId(roc_id)},{"$set": { "main_id": ObjectId(roc_id)}})
            except Exception, e:
                #raise Exception('Error when import ROC curve data: {}'.format(e))
                self.bind_object.logger.info("导入结果表%s时出错: %s" % (file, e))
            else:
                self.bind_object.logger.info("导入结果表%s成功!" % file)

    @report_check
    def add_roc_interval(self,roc_id,file):
        data_list=[]
        with open (file,"r") as f:
            lines=f.readlines()
            n = 1
            for line in lines[1:]:
                lst = line.strip('\r\n').split('\t')
                m=re.match(r'[0-9]*',lst[0])
                insert_data = {
                    'roc_id': ObjectId(roc_id),
                    'x': 1-float(m.group(0))/100.0,
                    'ymin': float(lst[1])/100.0,
                    'ymax': float(lst[3])/100.0,
                    "sort": n,
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
                n += 1
            try:
                collection = self.db["sg_roc_interval"]
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" % file)

    @report_check
    def add_roc_auc(self, roc_id, file, type):
        list_data = []
        dict = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lst = line.strip('\r\n').split('\t')
                dict['conf_level'] = float(lst[0])
                dict['min_auc'] = float(lst[1])
                dict['mean_auc'] = float(lst[2])
                dict['max_auc'] = float(lst[3])
            list_data.append(dict)
            try:
                if type in ['smooth']:
                    self.db['sg_roc'].update_one({'_id':ObjectId(roc_id)}, {'$set':{'auc_smooth':list_data}})
                else:
                    self.db['sg_roc'].update_one({"_id":ObjectId(roc_id)}, {'$set': {"auc": list_data}})
            except Exception, e:
                self.bind_object.logger.info("ROC主表导入auc出错：%s" % e)
            else:
                self.bind_object.logger.info("ROC主表导入auc成功！")

    @report_check
    def add_roc_best_loc(self, roc_id, file):
        best_loc = []
        tmplist = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lst = line.strip('\r\n').split('\t')
                tmplist["value"] =[(100-float(lst[2]))/100.0,(float(lst[5]))/100.0]
                tmplist["h_step_start"] = ((100-float(lst[2])) - (100-float(lst[3]))) / 100.0
                tmplist["h_step_end"] = ((100-float(lst[1])) - (100-float(lst[2]))) / 100.0
                tmplist["v_step_start"] = (float(lst[5]) - float(lst[4])) / 100.0
                tmplist["v_step_end"] = (float(lst[6]) - float(lst[5])) / 100.0
                tmplist["desc"] = "{}({},{})".format(float(lst[0]),(100-float(lst[2]))/100.0,(float(lst[5]))/100.0)
        best_loc.append(tmplist)
        try:
            self.db['sg_roc'].update_one({'_id': ObjectId(roc_id)}, {'$set': {'best_loc': best_loc,'status':'end'}})
        except Exception as e:
            self.bind_object.logger.info("ROC主表导入best_loc出错: %s" % e)
        else:
            self.bind_object.logger.info("ROC主表导入best_loc成功!")
