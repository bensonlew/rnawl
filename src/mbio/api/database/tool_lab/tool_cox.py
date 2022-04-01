# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20201012
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

class ToolCox(ApiBase):
    def __init__(self, bind_object):
        super(ToolCox, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_cox(self, main_id, cox, s3, params=None, project_sn='medical_transcriptome', task_id='medical_transcriptome'):
        # add main table info
        if main_id is None:
            name = "Cox" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Cox',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_tool_cox', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(cox, header=0, sep='\t')
        df['cox_id'] = main_id
        detail = df.to_dict('r')
        cox_columns = {'column': [{'field': 'feature', 'filter': False, 'sort': False, 'title': 'feature',
                                    'type': 'string'},
                                   {'field': 'beta', 'filter': False, 'sort': False,
                                    'title': 'beta',
                                    'type': 'float'},
                                   {'field': 'wald_test', 'filter': False, 'sort': False, 'title': 'wald_test',
                                    'type': 'float'},
                                   {'field': 'HR (95% CI for HR)', 'filter': False, 'sort': False, 'title': 'HR (95% CI for HR)',
                                    'type': 'float'},
                                   {'field': 'pvalue', 'filter': False, 'sort': False, 'title': 'pvalue',
                                    'type': 'float'}], 'condition': {}}
        column_data = json.dumps(cox_columns)
        self.create_db_table('sg_tool_cox_detail', detail)
        self.update_db_record('sg_tool_cox', main_id, status='end', main_id=main_id, s3_path=s3, column_data=column_data)
        return main_id


