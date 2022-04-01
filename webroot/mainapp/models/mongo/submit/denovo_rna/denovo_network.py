#!usr/bin/python
# -*- coding: utf-8 -*-
#__author__: konghualei
#last modify: 20161208
import datetime
from biocluster.config import Config
import json

class DenovoNetwork(object):
    def __init__(self):
        self.client = Config().mongo_client
        self.db_name = Config().MONGODB + '_rna'
        self.db = self.client[self.db_name]
    
    def add_network(self, params, softpower_png=None, module_png=None, name = None, module=None, project_sn = None, task_id = None):
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'gene_express_network_stat_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': "差异基因共表达网络主表",
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (json.dumps(params, sort_keys=True, separators=(',',':')) if isinstance(params, dict) else params),
            'module_png': module_png,
            'softpower_png': softpower_png,
            'status': 'start',
            'module': module
        }
        collection = self.db['sg_denovo_network']
        express_network_id = collection.insert_one(insert_data).inserted_id
        return express_network_id
   
