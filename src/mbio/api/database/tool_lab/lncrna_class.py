# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from collections import OrderedDict
import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import json
from api_base import ApiBase
import pandas as pd


class LncrnaClass(ApiBase):
    def __init__(self, bind_object):
        super(LncrnaClass, self).__init__(bind_object)

    @report_check
    def add_lncrna_class_detail(self, main_id, lncrna_class_detail):

        data_list = list()
        data_column = list()
        with open(lncrna_class_detail, 'r') as f:
            for line in f:
                cols = line.strip("\n").split("\t")
                if cols[1] != '0':
                    data_list.append({
                        'lncrna_class_id': main_id,
                        'lncrna_type': cols[0],
                        'num': int(cols[1])})
                    data_column.append({
                        'lncrna_class_id': main_id,
                        'name': cols[0],
                        'value': int(cols[1]),
                        'category': 'num',
                        'type': "column"
                    })

        table_dict = {
            "column": [
                {"field": "lncrna_type", "title": "Lncrna_Class", "filter": "false", "sort": "false", "type": "string"},
                {"field": "num", "title": "Num", "filter": "false", "sort": "false", "type": "string"}
            ],
            "condition": {}
        }
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))

        column_dict = {
            "name": "name",
            "data": "value",
            "condition": {"type": "column"}
        }
        column_info = json.dumps(column_dict, sort_keys=True, separators=(',', ':'))




        self.col_insert_data('lncrna_class_detail', data_list)
        self.col_insert_data('lncrna_class_column', data_column)


        self.update_db_record('lncrna_class',
                              query_dict={"main_id": ObjectId(main_id)},

                              update_dict={'status': 'end',
                                          'table_data': table_info,
                                          'column_data': column_info
                              })




if __name__ == '__main__':
    from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
    from biocluster.wsheet import Sheet
    import random


    data = {
        "id": "denovo_rna_v2_upgrade",
        #+ str(random.randint(1,10000)),
        #"id": "denovo_rna_v2",
        "project_sn": "denovo_rna_v2_upgrade",
        #+ str(random.randint(1,10000)),
        "type": "workflow",
        "name": "denovo_rna_v2.denovo_test_api",
        "options": {
        },
    }
    wsheet = Sheet(data=data)
    wf = DenovoTestApiWorkflow(wsheet)

    wf.IMPORT_REPORT_DATA = True
    wf.IMPORT_REPORT_AFTER_END = False
    wf.test_api = wf.api.api("denovo_rna_v2.cdslen")
