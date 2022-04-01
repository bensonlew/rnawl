# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

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
from api_base import ApiBase


class FeatureSelection(ApiBase):
    def __init__(self, bind_object):
        super(FeatureSelection, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_selection(self, main_id, result1=None, result2=None, params=None, project_sn='feature_selection', task_id='feature_selection'):
        # add main table info
        if main_id is None:
            name = "Feature_Selection" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Cox Feature Selection',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('feature_selection', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'status': 'end',
            'main_id': main_id,
            'pic_data': {
                'result1': result1,
                'result2': result2,
            },
        }
        self.update_db_record('feature_selection', query_dict=query_dict, update_dict=update_dict)
        return main_id
