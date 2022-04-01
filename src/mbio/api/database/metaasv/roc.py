# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os,re
import datetime
from bson.son import SON
import json
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check


class Roc(Base):
    def __init__(self, bind_object):
        super(Roc, self).__init__(bind_object)
        self._project_type = "metaasv"

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
            collection = self.db['roc']
            roc_id = collection.insert_one(insert_data).inserted_id
        except Exception, e:
            self.bind_object.logger.info("导入ROC主表出错: %s" % e)
        else:
            self.bind_object.logger.info("导入ROC主表成功!")
            return roc_id

    @report_check
    def add_roc_curve(self, roc_id, type, file):
        """
        导入线的数据为smooth和curve
        :param roc_id:
        :param type:
        :param file:
        :return:
        """
        data_list = []
        with open(file, "r") as f:
            lines = f.readlines()
            n = 1
            for line in lines[1:]:
                lst = line.strip().split("\t")
                insert_data = {
                    "x": (100-float(lst[0]))/100.0,
                    "y": float(lst[1])/100.0,
                    "type": type,
                    "name": "",
                    "roc_id": ObjectId(roc_id),
                    "sort": n,
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
                n += 1
            try:
                collection = self.db['roc_curve']
                collection.insert_many(data_list)
                main_collection = self.db["roc"]
                settled_params = {"software" : "R-3.3.1 (pROC)"}
                settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                if type in ["curve"]:
                    line_data = {
                        "line_data": {"name": "name",
                        "category": "",
                        "condition": {"type": "curve"}}
                        }
                    line_data_json = json.dumps(line_data, sort_keys=True, separators=(',', ':'))
                else:
                    line_data = {
                        "line_data": {"name": "name",
                        "category": "",
                        "condition": {"type": "smooth"}}
                        }
                    line_data_json = json.dumps(line_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": ObjectId(roc_id)},{"$set": { "main_id": ObjectId(roc_id),
                                                                            "settled_params": settled_params_json,
                                                                            "line_data": line_data_json,
                                                                            "auc": []
                                                                            }})
            except Exception, e:
                self.bind_object.logger.info("导入结果表%s时出错: %s" % (file, e))
            else:
                self.bind_object.logger.info("导入结果表%s成功!" % file)

    @report_check
    def add_roc_interval(self,roc_id,file):
        """
        导入区间图
        :param roc_id:
        :param file:
        :return:
        """
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
                    "name": "",
                    "type": "interval"
                }
                data_son = SON(insert_data)
                data_list.append(data_son)
                n += 1
            try:
                collection = self.db["roc_curve"]
                collection.insert_many(data_list)
                main_collection = self.db["roc"]
                line_data = {
                        "line_data": {"name": "name",
                        "category": "",
                        "condition": {"type": "interval"}}
                        }
                line_data_json = json.dumps(line_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id": ObjectId(roc_id)},{"$set": {"line_data": line_data_json,}})
            except Exception, e:
                self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
            else:
                self.bind_object.logger.info("导入%s结果表成功!" % file)

    @report_check
    def add_roc_auc(self, roc_id, file, type):
        """
        导入标签值（图上的标签）
        :param roc_id:
        :param file:
        :param type:
        :return:
        """
        dict = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lst = line.strip('\r\n').split('\t')
                dict['conf_level'] = float(lst[0])
                dict['min'] = float(lst[1])
                dict['mean'] = float(lst[2])
                dict['max'] = float(lst[3])
                dict['type'] = type
                dict['roc_id'] = ObjectId(roc_id)
        try:
            collection = self.db["roc_detail"]
            collection.insert_one(dict)
            main_collection = self.db["roc"]
            result = main_collection.find_one({"_id": ObjectId(roc_id)})
            auc = result["auc"]
            if type not in auc:
                auc.append(type)
            main_collection.update({"_id": ObjectId(roc_id)}, {"$set": {"auc": auc}})
        except Exception, e:
            self.bind_object.logger.info("roc_detail导入auc出错：%s" % e)
        else:
            self.bind_object.logger.info("roc_detail导入auc成功！")

    @report_check
    def add_roc_best_loc(self, roc_id, file, type):
        """
        导入临界值的最佳位置点
        :param roc_id:
        :param file:
        :param type:
        :return:
        """
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
                tmplist['roc_id'] = ObjectId(roc_id)
                tmplist['type'] = type
        best_loc.append(tmplist)
        try:
            collection = self.db["roc_detail"]
            collection.insert_many(best_loc)
        except Exception as e:
            self.bind_object.logger.info("roc_detail导入best_loc出错: %s" % e)
        else:
            self.bind_object.logger.info("roc_detail导入best_loc成功!")

    def add_roc_youden(self, roc_id, file, type):
        """
        导入约登指数
        :param roc_id: 主表id
        :param file: 文件路径
        :param type: 类型
        :return:
        """
        insert_data = {}
        with open(file, "r") as f:
            lines = f.readlines()
            line = lines[0].strip().split("\t")
            insert_data["value"] = line[0]
            insert_data['roc_id'] = ObjectId(roc_id)
            insert_data["type"] = type
        try:
            collection = self.db["roc_detail"]
            collection.insert_one(insert_data)
        except Exception as e:
            self.bind_object.logger.info("roc_detail导入youden_index出错: %s" % e)
        else:
            self.bind_object.logger.info("roc_detail导入youden_index成功!")
