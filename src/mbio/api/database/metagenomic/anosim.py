# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2017.11.06
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId


class Anosim(Base):
    def __init__(self, bind_object):
        super(Anosim, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_anosim(self, anosim_path, anosim_box_path, main=False, task_id=None, main_id=None,
                   anno_type=None, params=None, name=None):
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'Anosim_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'anno_type': anno_type
            }
            collection = self.db['anosim']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
        else:
            if main_id is None:
                self.bind_object.logger.error("main为False时需提供main_id!")
                self.bind_object.set_error("导表程序需要main_id参数", code="52800101")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        collection = self.db['anosim']
        with open(anosim_path + '/format_results.xls', 'r') as f:
            lines = f.readlines()
            line = lines[1].strip().split('\t')
        try:
            collection.update_one({"_id": ObjectId(main_id)},
                                  {"$set": {"statistic": round(float(line[1]),6), "p_value": round(float(line[2]),6), "permutation": int(line[3])}})
        except Exception as e:
            self.bind_object.logger.error("导入anosim%s信息出错:%s" % (anosim_path, e))
            self.bind_object.set_error("anosim导表出错", code="52800102")
        else:
            self.bind_object.logger.info("导入anosim%s信息成功!" % anosim_path)
        data_list = []
        with open(anosim_box_path + '/box_data.xls', 'r') as f:
            collection_detail = self.db["anosim_detail"]
            lines = f.readlines()
            for line in lines[1:]:
                data = [("anosim_id", main_id)]
                line = line.strip().split('\t')
                data.extend(
                    [("name", line[0]), ("min", line[1]), ("q1", line[2]), ("fliers", line[3]), ("median", line[4]),
                     ("q3", line[5]), ("max", line[6])])
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入anosim%s信息出错:%s" % (anosim_box_path, e))
            self.bind_object.set_error("anosim导表出错", code="52800102")
        else:
            self.bind_object.logger.info("导入anosim%s信息成功!" % anosim_box_path)
        return main_id
