# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import datetime
from biocluster.config import Config
import json


class DenovoMapping(object):
    def __init__(self):
        self.db_name = Config().MONGODB + '_rna'
        self.db = Config().mongo_client[self.db_name]

    def add_duplication_table(self, project_sn, task_id, name=None, params=None):
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "name": name if name else "duplication_origin",
            "status": "start",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_denovo_duplicate"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def add_rpkm_table(self, project_sn, task_id, name=None, params=None):
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "name": name if name else "rpkm_origin",
            "status": "start",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "curve_specimen": {"column1": "[0-0.3)", "column2": "[0.3-0.6)", "column3": "[0.6-3.5)", "column4": "[3.5-15)", "column5": "[15-60)", "column6": ">=60"},
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_denovo_rpkm"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def add_correlation_table(self, project_sn, task_id, name=None, params=None):
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "name": name if name else "correlation_origin",
            "status": "start",
            "desc": "",
            # "correlation_tree": correlation_tree,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_denovo_correlation"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id
