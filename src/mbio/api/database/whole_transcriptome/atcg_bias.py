# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import datetime
import json
import os
import unittest

import pandas as pd
from bson.objectid import ObjectId

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class AtcgBias(ApiBase):
    def __init__(self, bind_object):
        super(AtcgBias, self).__init__(bind_object)

    def creat_bias_main(self, collection_name, params=None, name=None, desc='atcg_bias_analyse'):
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn if self.bind_object else "small_rna",
            "task_id": self.bind_object.sheet.id if self.bind_object else "small_rna",
            "status": "end",
            "name": name if name else collection_name + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'desc': desc,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')) if params else None,
            'version': 'v1'
        }

        collection = self.db[collection_name]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    def export_first(self, result_dir, main_id, type_='known'):
        result_file = os.path.join(result_dir, '%s_first_bias_per.xls' % type_)
        result_pd = pd.read_table(result_file)
        result_pd.columns = ['length', 'a', 'g', 'c', 'u']
        result_pd = result_pd.fillna('_')
        result_pd['bias_id'] = ObjectId(main_id)
        result_pd['type'] = type_
        result_dict_list = result_pd.to_dict("records")
        self.create_db_table("bias_first", result_dict_list)

    def export_loc(self, result_dir, main_id, type_='known'):
        result_file = os.path.join(result_dir, '%s_loc_bias_per.xls' % type_)
        result_pd = pd.read_table(result_file)
        result_pd.columns = ['location', 'a', 'g', 'c', 'u']
        result_pd = result_pd.fillna('_')
        result_pd['bias_id'] = ObjectId(main_id)
        result_pd['type'] = type_
        result_dict_list = result_pd.to_dict("records")
        self.create_db_table("bias_location", result_dict_list)

    def run(self, family_path, params=None):
        main_id = self.creat_bias_main('bias', params=params)
        self.export_first(family_path, main_id, type_='known')
        self.export_first(family_path, main_id, type_='novel')
        self.export_first(family_path, main_id, type_='all')
        self.export_loc(family_path, main_id, type_='known')
        self.export_loc(family_path, main_id, type_='novel')
        self.export_loc(family_path, main_id, type_='all')
        self.db['bias'].update({"_id": main_id}, {"$set": {"main_id": main_id}})


class TestFunction(unittest.TestCase):
    def test_mongo(test):
        from mbio.workflows.prok_rna.prokrna_test_api import ProkrnaTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            "id": "dev_small_rna",
            # + str(random.randint(1,10000)),
            # "id": "denovo_rna_v2",
            "project_sn": "small_rna",
            # + str(random.randint(1,10000)),
            "type": "workflow",
            "name": "prok_rna.prokrna_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = ProkrnaTestApiWorkflow(wsheet)

        family_path = '/mnt/ilustre/users/sanger-dev/workspace/20181107/Single_AtcgBias15-58-31/AtcgBias/output'

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("small_rna.atcg_bias")
        params = {
            "name": 'atcg_analyse',
            "database": 'anminal',
        }

        wf.test_api.run(family_path, params=params)


if __name__ == '__main__':
    unittest.main()
