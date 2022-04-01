# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'fangyifei'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import types
import re
import datetime
import os
import json
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from bson.objectid import ObjectId


class FoldChange(ApiBase):
    def __init__(self, bind_object):
        super(FoldChange, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_fold_change_detail(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None, result_file=None):
        if main_id is None:
            name = "fold_change" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='fold_change',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('fold_change', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        data_list_bar = []
        with open(result_file) as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if float(item[2]) > 0:
                    category = "up"
                else:
                    category = "down"
                insert_data = {
                    "main_id": main_id,
                    "name": item[0],
                    "category": category,
                    "value": float(item[2]),
                    "type": "column",
                }
                data_list_bar.append(insert_data)
                data_list_bar = sorted(data_list_bar, key=lambda bar: bar["value"])
        self.create_db_table("fold_change_detail", data_list_bar)

        update_dict = json.dumps({
                "name": "name",
                "data": "value",
                "category": "category",
                "condition": {'type': "column"}
            })
        self.update_db_record("fold_change", main_id, main_id=main_id, column_data=update_dict, status="end")

