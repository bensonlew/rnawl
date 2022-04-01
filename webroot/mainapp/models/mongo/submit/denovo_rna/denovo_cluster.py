# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan'
# last_modify:20161011
import datetime
from biocluster.config import Config
import json


class DenovoExpress(object):
    def __init__(self):
        self.client = Config().mongo_client
        self.db_name = Config().MONGODB + '_rna'
        self.db = self.client[self.db_name]

    def add_cluster(self, params=None, name=None, project_sn=None, task_id=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '聚类分析主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(my_param, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'start'
        }
        collection = self.db['sg_denovo_cluster']
        cluster_id = collection.insert_one(insert_data).inserted_id
        return cluster_id
