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

class Duogene(ApiBase):
    def __init__(self, bind_object):
        super(Duogene, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_duogene(self, main_id, s3_path, params=None, project_sn='tool_lab', task_id='tool_lab',):
        # add main table info
        if main_id is None:
            name = "Duogene" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Duogene',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_duogene', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        self.update_db_record('sg_duogene', main_id, status='end', main_id=main_id, pdf_path = s3_path)



