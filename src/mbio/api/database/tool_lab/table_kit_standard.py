# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class TableKitStandard(ApiBase):
    def __init__(self, bind_object):
        super(TableKitStandard, self).__init__(bind_object)
        self._project_type = 'tool_lab'


    def add_standard(self, table_standard, main_id):
        # time_now = datetime.datetime.now()
        # name = 'Table_standard_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))

        # params = json.dumps({'task_id': task_id, 'submit_location': 'circrna', 'task_type': 2}, sort_keys=True)

        df = pd.read_table(table_standard)
        df_columns = list(df)
        # df_row = list(df[df_columns[0]])
        # df_col = df_columns[1:]
        df['standard_id'] = main_id
        # main_dict = {
        #     'task_id': task_id,
        #     'project_sn': project_sn,
        #     'name': name,
        #     'created_ts': time_now.strftime('%Y-%m-%d %H:%M:%S'),
        #     'desc': 'Table standard main table',
        #     'params': params,
        #     'status': 'start',
        # }
        # main_id = self.create_db_table('Table_standard', [main_dict])

        self.create_db_table('Table_standard_detail', df.to_dict('r'))
        self.update_db_record('Table_standard', main_id, insert_dict={'main_id': main_id, 'status': 'end', 'table_head':df_columns})
