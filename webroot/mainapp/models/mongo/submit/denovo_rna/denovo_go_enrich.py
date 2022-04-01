# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20161124
import datetime
from biocluster.config import Config
import json


class DenovoEnrich(object):
    def __init__(self):
        self.client = Config().mongo_client
        self.db_name = Config().MONGODB + '_rna'
        self.db = self.client[self.db_name]

    def add_go_enrich(self, name=None, params=None, project_sn=None, task_id=None, go_graph_dir=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'go_enrich_' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            'params': params,
            'status': 'end',
            'desc': 'go富集分析主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'go_directed_graph': '',
        }
        collection = self.db['sg_denovo_go_enrich']
        go_enrich_id = collection.insert_one(insert_data).inserted_id
        return go_enrich_id

    def add_go_regulate(self, name=None, params=None, project_sn=None, task_id=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'go_regulate_' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'start',
            'desc': 'go调控分析主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_go_regulate']
        go_regulate_id = collection.insert_one(insert_data).inserted_id
        return go_regulate_id

    def add_go_pval_sort(self, name=None, params=None, project_sn=None, task_id=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'GOEnrichRegulate_' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            'params': params,
            'status': 'end',
            'desc': 'go富集调控统计主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_go_pvalue']
        sort_id = collection.insert_one(insert_data).inserted_id
        return sort_id
