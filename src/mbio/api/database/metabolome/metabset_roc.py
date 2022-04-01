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
import pandas as pd
#from biocluster.config import Config

class MetabsetRoc(Base):
    def __init__(self, bind_object):
        super(MetabsetRoc, self).__init__(bind_object)
        self._project_type = "metabolome"
        self.metabset_roc_table_data = dict()

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
            collection = self.db['metabset_roc']
            roc_id = collection.insert_one(insert_data).inserted_id
        except Exception as e:
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
            ori_metab= ''
            for line in lines[1:]:
                lst = line.strip("\r\n").split("\t")
                if lst[0] != ori_metab:
                    n = 1
                    ori_metab = lst[0]
                insert_data = {
                    "x": (100-float(lst[1]))/100.0,
                    "y": float(lst[2])/100.0,
                    "type": type,
                    "roc_id": ObjectId(roc_id),
                    "sort": n,
                    "metab_id" : lst[0]
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
                n += 1
            try:
                collection = self.db['metabset_roc_curve']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.info("导入结果表%s时出错: %s" % (file, e))
            else:
                self.bind_object.logger.info("导入结果表%s成功!" % file)

    @report_check
    def add_roc_interval(self,roc_id,file):
        data_list=[]
        with open (file,"r") as f:
            lines=f.readlines()
            n = 1
            ori_metab = ''
            for line in lines[1:]:
                lst = line.strip('\r\n').split('\t')
                m=re.match(r'([0-9]*)',lst[1])
                if lst[0] != ori_metab:
                    n = 1
                    ori_metab = lst[0]
                insert_data = {
                    'roc_id': ObjectId(roc_id),
                    'x': 1-float(m.group(1))/100.0,
                    'ymin': float(lst[2])/100.0,
                    'ymax': float(lst[4])/100.0,
                    "sort": n,
                    "metab_id" : lst[0]   #add v3 20200325
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
                n += 1
            try:
                collection = self.db["metabset_roc_interval"]
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" % file)

    def change_float_100(self,instr):
        try:
            return round(float(instr)/100,4)
        except:
            return instr

    @report_check
    def add_roc_auc(self, roc_id, file, type):

        if type in ['smooth']:
            mk = 'auc_smooth'
        else:
            mk = 'auc'
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                dictt = {}
                lst = line.strip('\r\n').split('\t')
                dictt['conf_level'] = self.change_float_100(lst[1])
                dictt['min_auc'] = self.change_float_100(lst[2])
                dictt['mean_auc'] = self.change_float_100(lst[3])
                dictt['max_auc'] = self.change_float_100(lst[4])

                if lst[0] not in self.metabset_roc_table_data:
                    self.metabset_roc_table_data[lst[0]] = {'roc_id' : ObjectId(roc_id),'metab_id': lst[0]}

                self.metabset_roc_table_data[lst[0]][mk] =  dictt
        self.bind_object.logger.info("%s file 处理完成"%file)

    @report_check
    def add_roc_best_loc(self, roc_id, file):
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                tmplist = {}
                lst = line.strip('\r\n').split('\t')
                # tmplist["value"] =[(100-float(lst[2]))/100.0,(float(lst[5]))/100.0]
                # tmplist["h_step_start"] = ((100-float(lst[2])) - (100-float(lst[3]))) / 100.0
                # tmplist["h_step_end"] = ((100-float(lst[1])) - (100-float(lst[2]))) / 100.0
                # tmplist["v_step_start"] = (float(lst[5]) - float(lst[4])) / 100.0
                # tmplist["v_step_end"] = (float(lst[6]) - float(lst[5])) / 100.0
                # tmplist["desc"] = "{}({},{})".format(float(lst[0]),(100-float(lst[2]))/100.0,(float(lst[5]))/100.0)
                tmplist["value"] =[(100-float(lst[3]))/100.0,(float(lst[6]))/100.0]
                tmplist["h_step_start"] = ((100-float(lst[3])) - (100-float(lst[4]))) / 100.0
                tmplist["h_step_end"] = ((100-float(lst[2])) - (100-float(lst[3]))) / 100.0
                tmplist["v_step_start"] = (float(lst[6]) - float(lst[5])) / 100.0
                tmplist["v_step_end"] = (float(lst[7]) - float(lst[6])) / 100.0
                tmplist["desc"] = "{}({},{})".format(float(lst[1]),(100-float(lst[3]))/100.0,(float(lst[6]))/100.0)

                if lst[0] not in self.metabset_roc_table_data:
                    self.metabset_roc_table_data[lst[0]] = {'roc_id' : ObjectId(roc_id),'metab_id': lst[0]}

                self.metabset_roc_table_data[lst[0]]['best_loc'] = tmplist

        self.bind_object.logger.info("%s file 处理完成"%file)

    @report_check
    def add_roc_table(self,roc_id, metab_desc):
        self.bind_object.logger.info('运行add_roc_auc和add_roc_best_loc成功，然后运行add_roc_table')
        desc = pd.read_table(metab_desc, sep='\t',index_col=0)
        add_info = {
            'Metabolite' : 'metab_name'
        }

        if 'ID' in desc.columns:
            var_head = 'o_id'
            add_info['ID'] = 'o_id'
        else:
            var_head = ''

        for metab_id in self.metabset_roc_table_data:
                for ori_k in add_info:
                    mongo_k = add_info[ori_k]
                    info = desc.loc[metab_id][ori_k]
                    self.metabset_roc_table_data[metab_id][mongo_k] = info

        self.db['metabset_roc_table'].insert_many(self.metabset_roc_table_data.values())
        if var_head:
            self.db['metabset_roc'].update_one({'_id': ObjectId(roc_id)}, {'$set': {'status':'end', 'var_head': var_head}})
        else:
            self.db['metabset_roc'].update_one({'_id': ObjectId(roc_id)}, {'$set': {'status':'end'}})
        self.bind_object.logger.info('导入add_roc_table成功')

