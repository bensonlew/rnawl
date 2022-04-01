# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20210202
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class Busco(ApiBase):
    def __init__(self, bind_object):
        super(Busco, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_busco(self, main_id, s3_path, params=None, project_sn='tool_lab', task_id='tool_lab',):
        # add main table info
        if main_id is None:
            name = "Busco" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Busco',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_busco', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        # df = pd.read_table(summary_busco, header=0, sep='\t')
        # columns_list = list()
        # columns_list.append({'field': 'name', 'filter': False, 'sort': False, 'title': 'Name', 'type': 'string'})
        # columns_list.append({'field': 'number', 'filter': False, 'sort': False, 'title': 'Number', 'type': 'int'})
        # columns_list.append({'field': 'persent', 'filter': False, 'sort': False, 'title': 'Persent', 'type': 'float'})
        # data_columns = {'column': columns_list, 'condition': {}}
        # columns_data = json.dumps(data_columns)
        # df['busco_id'] = main_id
        # detail = df.to_dict('r')
        # self.create_db_table('sg_busco_detail', detail)
        # self.update_db_record('sg_busco', main_id, main_id=main_id, column_data_detail=columns_data)
        # df = pd.read_table(summary_busco, header=0, sep='\t', usecols=['name', 'persent'])
        # df = df[1:-1]
        # df.rename(columns={'name': 'category', 'persent': 'value'}, inplace=True)
        # df['name'] = sample_name
        # df["type"] = "column"
        # df['busco_id'] = main_id
        # stacked_detail = df.to_dict('r')
        # self.create_db_table('sg_busco_stacked_detail', stacked_detail)
        # dict_a = {"name": "name", "data": "value", "category": "category", "condition": {"type": "column"}}
        # dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        # self.update_db_record('sg_busco', main_id, main_id=main_id, column_data_stacked=dict_b)
        #
        # text_detail = [{'text': text, 'busco_id': main_id, 'name': sample_name}]
        # self.create_db_table('sg_busco_text_detail', text_detail)
        # text_dict = {"name": "name", "condition": {"type": "text"}}
        # text_json = json.dumps(text_dict, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_busco', main_id, status='end', main_id=main_id, pdf_path = s3_path)



