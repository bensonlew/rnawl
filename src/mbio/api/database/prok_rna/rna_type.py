# !/usr/bin/python
# -*- coding: utf-8 -*-
import types
import os
import json
import unittest
import datetime
import sqlite3
import re
from collections import OrderedDict, defaultdict
from bson.son import SON
from bson.objectid import ObjectId
from mbio.api.database.prok_rna.api_base import ApiBase
from biocluster.api.database.base import Base, report_check
import glob
import copy


class RnaType(ApiBase):
    def __init__(self, bind_object):
        super(RnaType, self).__init__(bind_object)
        self._project_type = 'prok_rna'

    @report_check
    def add_rna_type(self, matrix, params=None, name=None, desc='rna_type'):
        """
        添加rna_type的主表函数
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn if self.bind_object else "prok_rna",
            "task_id": self.bind_object.sheet.id if self.bind_object else "prok_rna",
            "status": "end",
            "name": name if name else "Rna_type_" + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'desc': desc,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')) if params else "rna_type"
        }

        collection = self.db['sg_rna_type']
        inserted_id = collection.insert_one(insert_data).inserted_id
        self.add_rnatype_table_detail(matrix, inserted_id)

    @report_check
    def add_rnatype_table_detail(self, matrix, rna_type_id):
        """
        添加rna_type的详情表函数
        """
        with open(matrix, "r") as f:
            for line in f:
                items = line.strip().split("\t")
                insert_data = {
                    'rna_type_id': rna_type_id,
                    'seq_id': items[0],
                    'type': items[1]
                }
                collection = self.db['sg_rna_type_detail']
                collection.insert_one(insert_data)


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """
    def test_mongo(test):
        from mbio.workflows.prok_rna.prokrna_test_api import ProkrnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            "id": "rna_type",
            "project_sn": "prok_rna_srna",
            "type": "workflow",
            "name": "prok_rna.prokrna_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = ProkrnaTestApiWorkflow(wsheet)
        matrix = '/mnt/ilustre/users/sanger-dev/workspace/20190402/Prokrna_tsg_33758/ExtractBiotype/output/gene_biotype.txt'

        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("prok_rna.rna_type")
        wf.test_api.add_rna_type(matrix)

if __name__ == '__main__':
    unittest.main()