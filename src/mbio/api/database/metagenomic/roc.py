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
from bson.objectid import ObjectId

class Roc(Base):
    def __init__(self, bind_object):
        super(Roc, self).__init__(bind_object)
        self._project_type = "metagenomic"

    @report_check
    def add_roc(self,params=None,name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': 'ROC 分析',
            'created_ts': created_ts,
            'name': name if name else "ROC",
            'params': params,
            'status': 'end'
        }
        try:
            collection = self.db['roc']
            roc_id = collection.insert_one(insert_data).inserted_id
        except Exception, e:
            #raise Exception('导入roc主表异常:{}'.format(e))
            self.bind_object.set_error('导入roc主表异常:%s', variables=(e), code="52803201")
        return roc_id

    @report_check
    def add_roc_curve(self, roc_id, type,file):
        data_list=[]
        with open (file,"r") as f:
            lines=f.readlines()
            n = 1
            for line in lines[1:]:
                lin = line.strip('\r\n').split('\t')
                insert_data = {
                    'roc_id':ObjectId(roc_id),
                    'type':type,
                    'x': (100-float(lin[0]))/100.0,
                    'y': float(lin[1])/100.0,
                    "sort":n,
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
                n +=1
            try:
                collection = self.db["roc_curve"]
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" % file)

    @report_check
    def add_roc_interval(self,roc_id,file):
        data_list=[]
        with open (file,"r") as f:
            lines=f.readlines()
            n=1
            for line in lines[1:]:
                lin = line.strip('\r\n').split('\t')
                m=re.match(r'[0-9]*',lin[0])
                insert_data = {
                    'roc_id':ObjectId(roc_id),
                    'x': 1-float(m.group(0))/100.0,
                    'ymin': float(lin[1])/100.0,
                    'ymax': float(lin[3])/100.0,
                    'sort':n,
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
                n +=1
            try:
                collection = self.db["roc_interval"]
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" % file)

    @report_check
    def add_roc_auc(self,roc_id,file,type):
        list = []
        dict={}
        with open (file,"r") as f:
            lines=f.readlines()
            for line in lines[1:]:
                lin = line.strip('\r\n').split('\t')
                dict['conf_level'] =lin[0]
                dict['min_auc'] = lin[1]
                dict['mean_auc'] = lin[2]
                dict['max_auc'] = lin[3]
        list.append(dict)
        try:
            if type in ['smooth']:
                self.db['roc'].update_one({'_id': ObjectId(roc_id)}, {'$set': {'auc_smooth': list}})
            else:
                self.db['roc'].update_one({'_id': ObjectId(roc_id)}, {'$set': {'auc': list}})
        except Exception as e:
            #raise Exception('更新roc主表file出错:{}'.format(e))
            self.bind_object.set_error('更新roc主表file出错:%s', variables=(e), code="52803202")

    @report_check
    def add_roc_best_loc(self, roc_id, file):
        best_loc = []
        tmplist = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip('\r\n').split('\t')
                tmplist["value"] =[(100-float(lin[2]))/100.0,(float(lin[5]))/100.0]
                tmplist["h_step_start"] = ((100-float(lin[2])) - (100-float(lin[3]))) / 100.0
                tmplist["h_step_end"] = ((100-float(lin[1])) - (100-float(lin[2]))) / 100.0
                tmplist["v_step_start"] = (float(lin[5]) - float(lin[4])) / 100.0
                tmplist["v_step_end"] = (float(lin[6]) - float(lin[5])) / 100.0
                tmplist["desc"] = "{}({},{})".format(float(lin[0]),(100-float(lin[2]))/100.0,(float(lin[5]))/100.0)
        best_loc.append(tmplist)
        try:
            self.db['roc'].update_one({'_id': ObjectId(roc_id)}, {'$set': {'best_loc': best_loc,'status':'end'}})
        except Exception as e:
            #raise Exception('更新roc主表file出错:{}'.format(e))
            self.bind_object.set_error('更新roc主表file出错:%s', variables=(e), code="52803203")




