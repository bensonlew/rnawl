# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Km(ApiBase):
    def __init__(self, bind_object):
        super(Km, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_km(self, project_sn='tool_lab', main_id=None, task_id='tool_lab', s3_output=None):
        if main_id is None:
            name = "Km_survival" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='km main table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('km_survival', [main_info])

        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        self.update_db_record('km_survival', main_id, status="end", main_id=main_id, pic_data=s3_output)
        return main_id
