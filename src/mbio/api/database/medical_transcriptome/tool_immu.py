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
from mbio.api.database.medical_transcriptome.api_base import ApiBase

class ToolImmu(ApiBase):
    def __init__(self, bind_object):
        super(ToolImmu, self).__init__(bind_object)

    def add_immu(self, main_id, immu, method,params=None, project_sn='medical_transcriptome', task_id='medical_transcriptome'):
        # add main table info
        if main_id is None:
            name = "Immunedeconv" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Immunedeconv',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_tool_immunedeconv', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(immu, header=0, sep='\t')
        row_name = df['cell_type'].tolist()
        if method.lower() in ['xcell', 'cibersort']:
            cell_type = row_name[:-3]
        else:
            cell_type = row_name
        df['immu_id'] = main_id
        detail = df.to_dict('r')
        self.create_db_table('sg_tool_immunedeconv_detail', detail)
        self.update_db_record('sg_tool_immunedeconv', main_id, status='end', main_id=main_id, cell_type=cell_type, row_name=row_name)
        return main_id


