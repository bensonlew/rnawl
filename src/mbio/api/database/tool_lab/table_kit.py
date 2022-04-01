# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class TableKit(ApiBase):
    def __init__(self, bind_object):
        super(TableKit, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_table_standard(self, project_sn='tool_lab', main_id=None, task_id='tool_lab'):
        if main_id is None:
            name = "Standard" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")

            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Standard main table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('Table_standard', [main_info])

        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        self.update_db_record('Table_standard', main_id, status="end", main_id=main_id)
        return main_id
