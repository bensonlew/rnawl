# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import datetime
import json
import os
import unittest

import pandas as pd
from bson.objectid import ObjectId

from mbio.api.database.whole_transcriptome.api_base import ApiBase


class SmallrnaFamilyAnalyse(ApiBase):
    def __init__(self, bind_object):
        super(SmallrnaFamilyAnalyse, self).__init__(bind_object)

    def creat_family_main(self, collection_name, params=None, name=None, desc='family_analyse'):
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

    def export_family(self, result_dir, main_id):
        result_file = os.path.join(result_dir, 'family.species.xls')
        result_pd = pd.read_table(result_file, index_col=False)
        result_pd.rename(columns={result_pd.columns[0]: "fullname", result_pd.columns[1]: 'abbreviation'}, inplace=True)
        result_pd = result_pd.fillna('_')
        result_pd['family_id'] = ObjectId(main_id)
        result_dict_list = result_pd.to_dict("records")
        # for i in result_dict_list:
        #     print(i['mir-96'])
        self.create_db_table("family_stat", result_dict_list)

    def export_known(self, result_dir, main_id):
        result_file = os.path.join(result_dir, 'known_miR_family.xls')
        result_pd = pd.read_table(result_file)
        result_pd = result_pd.fillna('_')
        result_pd.columns = ['mirna', 'mir_pre', 'family']
        result_pd['family_id'] = ObjectId(main_id)
        result_pd['type'] = 'known'
        result_dict_list = result_pd.to_dict("records")
        for r_d in result_dict_list:
            r_d['mir_pre_list'] = r_d['mir_pre'].split(';')
            r_d['family_list'] = r_d['family'].split(';')
        self.create_db_table("family_detail", result_dict_list)

    def export_novel(self, result_dir, main_id):
        result_file = os.path.join(result_dir, 'novel_miR_family.xls')
        result_pd = pd.read_table(result_file)
        result_pd = result_pd.fillna('_')
        result_pd.columns = ['mirna', 'family']
        result_pd['mir_pre'] = '-'
        result_pd['family_id'] = ObjectId(main_id)
        result_pd['type'] = 'novel'
        result_dict_list = result_pd.to_dict("records")
        for r_d in result_dict_list:
            r_d['mir_pre_list'] = r_d['mir_pre'].split(';')
            r_d['family_list'] = r_d['family'].split(';')
        self.create_db_table("family_detail", result_dict_list)

    def run(self, family_path, params=None):
        main_id = self.creat_family_main('family', params=params)
        self.export_family(family_path, main_id)
        self.export_known(family_path, main_id)
        self.export_novel(family_path, main_id)
        with open(os.path.join(family_path, 'family.species.xls'), 'r') as fam_r:
            family_names = ';'.join(fam_r.readline().strip().split('\t')[2:])
        self.db['family'].update({"_id": main_id}, {"$set": {"main_id": main_id, "family_names": family_names}})


class TestFunction(unittest.TestCase):
    def test_mongo(test):
        from mbio.workflows.prok_rna.prokrna_test_api import ProkrnaTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            "id": "small_rna1",
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

        family_path = '/mnt/ilustre/users/sanger-dev/workspace/20181220/Smallrna_t33076/SmallrnaFamilyAnalyse/output'

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("small_rna.smallrna_family_analyse")
        params = {
            "name": 'family_analyse',
            "database": 'anminal',
        }

        wf.test_api.run(family_path, params=params)


if __name__ == '__main__':
    unittest.main()
