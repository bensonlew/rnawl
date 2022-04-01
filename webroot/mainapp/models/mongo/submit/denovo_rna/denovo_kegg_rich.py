# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify:20161124
import datetime
from biocluster.config import Config
import json
from mainapp.models.mongo.core.base import Base


class DenovoKeggRich(Base):
    def __init__(self):
        super(DenovoKeggRich, self).__init__(bind_object)
        self._project_type = 'ref_rna'


    def add_kegg_rich(self, name=None, params=None, project_sn=None, task_id=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'KeggEnrich_' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'desc': 'kegg富集分析',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_kegg_enrich']
        enrich_id = collection.insert_one(insert_data).inserted_id
        return enrich_id

    def add_kegg_regulate(self, name=None, params=None, project_sn=None, task_id=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'KeggRegulate_' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'desc': 'kegg调控统计分析',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_kegg_regulate']
        regulate_id = collection.insert_one(insert_data).inserted_id
        return regulate_id

    def add_kegg_pval_sort(self, name=None, params=None, project_sn=None, task_id=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'KeggEnrichRegulate_' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'status': 'end',
            'desc': 'kegg富集调控pval排序',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db['sg_denovo_kegg_pvalue']
        sort_id = collection.insert_one(insert_data).inserted_id
        return sort_id
