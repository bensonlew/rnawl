# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20190902

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
import pandas as pd
import numpy as np


class PanVariation(Base):
    def __init__(self, bind_object):
        super(PanVariation, self).__init__(bind_object)
        self._project_type = "bac_comparative"
        #self.id = 'tsg_123'
        #self.project_sn = '188_5b5acb3018'

    @report_check
    def add_pan(self, params=None, name=None):
        """
        pan_genomes的主表；params中记录了group_detail等字段
        :param params:
        :return:
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "同源基因突变分析主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Pan_var_Origin",
        }
        collection = self.db["pan_var"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}})
        return inserted_id

    def add_pan_variation_stat(self, inserted_id, variation):
        """
        导入变异分析统计表
        :param inserted_id: 主表ID
        :param variation: 导入变异分析统计的表格
        :return:
        """
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        with open(variation, 'r') as f:
            f.readline()
            new_insert = 0
            for log_index,line in enumerate(f):
                if log_index % 200000 == 0 and log_index > 0:
                    new_insert = 1
                    self.bind_object.logger.info("have done %s lines" % log_index)
                line = line.strip().split('\t')
                dn = int(line[2])
                ds = int(line[3])
                if ds == 0 or dn ==0:
                    stat = 0
                else:
                    stat = float(dn) / ds
                data = {
                    "cluster_id": line[0],
                    "var_id": new_inserted_id,
                    "indel_num": int(line[1]),
                    "nonsy_num": int(line[2]),
                    "sy_num": int(line[3]),
                    "dn_ds": stat
                    }
                data_list.append(data)
                if new_insert == 1 and log_index > 0:
                    self.bind_object.logger.info("开始导入数据库前%s lines" % log_index)
                    try:
                        collection = self.db['pan_variation']
                        collection.insert_many(data_list,ordered=True)
                        new_insert = 0
                        data_list = []
                    except Exception,e:
                        self.bind_object.logger.error("导入前%slines失败:%s" % (log_index,e))
                        self.bind_object.set_error("导入数据失败")
            self.bind_object.logger.info("导入剩余数据")
            try:
                collection = self.db["pan_var_detail"]
                collection.insert_many(data_list)
                main_collection = self.db["pan_var"]
                main_collection.update({'_id': new_inserted_id},{'$set':{'main_id': new_inserted_id}})
            except Exception, e:
                self.bind_object.logger.info("导入variation结果表出错:%s" % (e))
            else:
                self.bind_object.logger.info("导入variation结果表成功!")


