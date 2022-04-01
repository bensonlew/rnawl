# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last_modify:20160919
import datetime
from biocluster.config import Config
import json


class DenovoExpress(object):
    def __init__(self):
        self.client = Config().mongo_client
        self.db_name = Config().MONGODB + '_rna'
        self.db = self.client[self.db_name]

    def add_express(self, rsem_dir=None, samples=None, params=None, name=None, bam_path=None, project_sn=None, task_id=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'express_matrix_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': '表达量计算主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'specimen': samples,
            'status': 'start',
            'bam_path': bam_path,
        }
        collection = self.db['sg_denovo_express']
        express_id = collection.insert_one(insert_data).inserted_id
        return express_id

    def add_express_diff(self, params, samples=None, compare_column=None, name=None, project_sn=None, task_id=None, samples_detail=None, express_id=None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'gene_express_diff_stat_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': '表达量差异检测主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'specimen': samples,
            'status': 'start',
            'compare_column': compare_column,
            'group_detail': samples_detail,
            'express_id': express_id
        }
        collection = self.db['sg_denovo_express_diff']
        express_diff_id = collection.insert_one(insert_data).inserted_id
        return express_diff_id
