# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last modify: 2018.04.2

from biocluster.api.database.base import Base, report_check
import datetime
import json
from bson.son import SON
import pandas as pd
from collections import Counter
from bson.objectid import ObjectId


class CircosTable(Base):
    def __init__(self, bind_object):
        super(CircosTable, self).__init__(bind_object)
        self._project_type = "bacgenome"

    def add_circos_table(self, params=None, project_sn=None, task_id=None):
        task_id = task_id if task_id else self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "circos可选分析项目",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
        }
        collection = self.db["circos_table"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set':{'main_id':inserted_id}}) #guanqing 20180813
        self.update_table(task_id=task_id, main_id=inserted_id)
        return inserted_id

    @report_check
    def update_table(self, task_id=None, main_id=None):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        detail_table = self.db['circos_table_detail']
        task_id = task_id if task_id else self.bind_object.sheet.id
        try:
            trna = self.db['trna_predict']
            trna_id = trna.find_one({"task_id": task_id})['_id']
            trna_detail = self.db['trna_predict_detail']
            for one in trna_detail.find({"predict_id": trna_id}):
                if one['location'].startswith('Scaffold'):
                    detail_table.update_one(
                        {"circos_table_id": main_id, "specimen_id": one['specimen_id']},
                        {'$set': {'ncRNA': 1}}, upsert=True)
                else:
                    detail_table.update_one(
                        {"circos_table_id": main_id, "location": one['location'], "specimen_id": one['specimen_id']},
                        {'$set': {'ncRNA': 1}}, upsert=True)
        except Exception as e:
            self.bind_object.logger.info('没有trna信息')
        else:
            self.bind_object.logger.info('导入trna信息成功')
        try:
            rrna = self.db['rrna_predict']
            rrna_id = rrna.find_one({"task_id": task_id})['_id']
            rrna_detail = self.db['rrna_predict_detail']
            for one in rrna_detail.find({"predict_id": rrna_id}):
                if one['location'].startswith('Scaffold'):
                    detail_table.update_one(
                        {"circos_table_id": main_id, "specimen_id": one['specimen_id']},
                        {'$set': {'ncRNA': 1}}, upsert=True)
                else:
                    detail_table.update_one(
                        {"circos_table_id": main_id, "location": one['location'], "specimen_id": one['specimen_id']},
                        {'$set': {'ncRNA': 1}}, upsert=True)
        except Exception as e:
            self.bind_object.logger.info('没有rrna信息')
        else:
            self.bind_object.logger.info('导入rrna信息成功')
        try:
            type_main = self.db['prephage']
            type_main_id = type_main.find_one({"task_id": task_id})['_id']
            type_main_detail = self.db['prephage_detail']
            for one in type_main_detail.find({"prephage_id": type_main_id}):
                if one['location'].startswith('Scaffold'):
                    detail_table.update_one(
                        {"circos_table_id": main_id, "specimen_id": one['specimen_id']},
                        {'$set': {'Prephage': 1}}, upsert=True)
                else:
                    detail_table.update_one(
                        {"circos_table_id": main_id, "location": one['location'], "specimen_id": one['specimen_id']},
                        {'$set': {'Prephage': 1}}, upsert=True)
        except Exception as e:
            self.bind_object.logger.info('没有prephage信息')
        else:
            self.bind_object.logger.info('导入prephage信息成功')
        try:
            type_main = self.db['island']
            type_main_id = type_main.find_one({"task_id": task_id})['_id']
            type_main_detail = self.db['island_detail']
            for one in type_main_detail.find({"island_id": type_main_id}):
                if one['location'].startswith('Scaffold'):
                    detail_table.update_one(
                        {"circos_table_id": main_id, "specimen_id": one['specimen_id']},
                        {'$set': {'GI': 1}}, upsert=True)
                else:
                    detail_table.update_one(
                        {"circos_table_id": main_id, "location": one['location'], "specimen_id": one['specimen_id']},
                        {'$set': {'GI': 1}}, upsert=True)
        except Exception as e:
            self.bind_object.logger.info('没有基因组岛信息')
        else:
            self.bind_object.logger.info('导入基因组岛信息成功')
        try:
            type_main = self.db['is']
            type_main_id = type_main.find_one({"task_id": task_id})['_id']
            type_main_detail = self.db['is_detail']
            for one in type_main_detail.find({"is_id": type_main_id}):
                if one['location'].startswith('Scaffold'):
                    detail_table.update_one(
                        {"circos_table_id": main_id, "specimen_id": one['sample']},
                        {'$set': {'IS': 1}}, upsert=True)
                else:
                    detail_table.update_one(
                        {"circos_table_id": main_id, "location": one['location'], "specimen_id": one['sample']},
                        {'$set': {'IS': 1}}, upsert=True)
        except Exception as e:
            self.bind_object.logger.info('没有插入序列信息')
        else:
            self.bind_object.logger.info('导入插入序列信息')
        try:
            type_main = self.db['integron']
            type_main_id = type_main.find_one({"task_id": task_id})['_id']
            type_main_detail = self.db['integron_detail']
            for one in type_main_detail.find({"integron_id": type_main_id}):
                if one['location'].startswith('Scaffold'):
                    detail_table.update_one(
                        {"circos_table_id": main_id, "specimen_id": one['sample']},
                        {'$set': {'Integron': 1}}, upsert=True)
                else:
                    detail_table.update_one(
                        {"circos_table_id": main_id, "location": one['location'], "specimen_id": one['sample']},
                        {'$set': {'Integron': 1}}, upsert=True)
        except Exception as e:
            self.bind_object.logger.info('没有整合子信息')
        else:
            self.bind_object.logger.info('导入整合子信息成功')
