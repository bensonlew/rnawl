# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

import logging
import sys
import time
import json
from biocluster.config import Config
from concurrent.futures import ThreadPoolExecutor

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


class AddbatchMongo(object):
    def __init__(self, project_type):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self._project_type = project_type

    def run(self):
        self.add_params_batch()

    def add_params_batch(self):
        sg_diff = self.db['sg_diff']
        document = sg_diff.find({})
        for record in document:
            try:
                record['is_batch']
                continue
            except:
                try:
                    params = record['params']
                    params = json.loads(params)
                    params_old = params
                    params['is_batch'] = "False"
                    params_new = json.dumps(params, sort_keys=True, separators=(',', ':'))
                    params_old = json.dumps(params_old, sort_keys=True, separators=(',', ':'))
                    sg_diff.update({'_id': record["_id"]}, {'$set': {'params': params_new}})
                except:
                    continue

if __name__ == '__main__':
    logging.info('Usage: python add_batch.py <project_type>')
    # task_id = sys.argv[1]
    project_type = sys.argv[1]
    start_time = time.time()
    inst = AddbatchMongo(project_type)
    inst.run()
    end_time = time.time()
    logging.info('Elapsed time: {}'.format(end_time - start_time))
