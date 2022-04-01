# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import datetime
from biocluster.config import Config
import json


class GeneStructure(object):
    def __init__(self):
        self.db_name = Config().MONGODB + '_rna'
        self.db = Config().mongo_client[self.db_name]

    def add_ssr_table(self, project_sn, task_id, name=None, params=None):
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "name": name if name else "ssr_origin",
            "status": "start",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_denovo_ssr"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id
