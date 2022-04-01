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

class Translation(ApiBase):
    def __init__(self, bind_object):
        super(Translation, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_translation(self, main_id, translation_file, params=None, project_sn='tool_lab', task_id='tool_lab',):
        # add main table info
        if main_id is None:
            name = "Plant_translation" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Plant Name Translation',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_plant_trans', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(translation_file, header=0, sep='\t')
        columns_list = list()
        columns_list.append({'field': 'Chinese name', 'filter': False, 'sort': False, 'title': '中文名称', 'type': 'string'})
        columns_list.append({'field': 'English name', 'filter': False, 'sort': False, 'title': '英文名称', 'type': 'string'})
        columns_list.append({'field': 'Latin name', 'filter': False, 'sort': False, 'title': '拉丁学名', 'type': 'string'})
        data_columns = {'column': columns_list, 'condition': {}}
        columns_data = json.dumps(data_columns)
        df['translation_id'] = main_id
        detail = df.to_dict('r')
        if detail:
            self.create_db_table('sg_plant_trans_detail', detail)
        self.update_db_record('sg_plant_trans', main_id, main_id=main_id, column_data_detail=columns_data, status='end')



