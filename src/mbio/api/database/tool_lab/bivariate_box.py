# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'fangyifei'

import os
from bson.objectid import ObjectId
import pandas as pd
import glob
from mbio.api.database.whole_transcriptome.api_base import ApiBase
import json
import numpy as np
import copy
import datetime


class BivariateBox(ApiBase):
    def __init__(self, bind_object):
        super(BivariateBox, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_bivariate_box_detail(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None, box_file=None, scatter_file=None):
        if main_id is None:
            name = "bivariate_box" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='bivariate_box',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('bivariate_box', [main_info])
        else:
            main_id = ObjectId(main_id)
            # self.update_db_record('bivariate_box')
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        data_list1 = []
        with open(box_file) as f1:
            lines = f1.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                if item[2] == "":
                    item[2] = item[3]
                if len(item) == 6:
                    item.append(item[5])
                elif item[6] == "":
                    item[6] = item[5]
                insert_data1 = {
                    "main_id": main_id,
                    "name": item[0],
                    "box_name": item[1],
                    "category": item[1],
                    "min": float(item[2]),
                    "q1": float(item[3]),
                    "median": float(item[4]),
                    "q3": float(item[5]),
                    "max": float(item[6]),
                    "type": "box"
                }
                data_list1.append(insert_data1)
        self.create_db_table('bivariate_box_detail', data_list1)

        # 离群值
        data_list2 = []
        with open(scatter_file) as f2:
            lines = f2.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data2 = {
                    "main_id": main_id,
                    "sample": item[0],
                    "x": item[1],
                    "y": float(item[3]),
                    "category": item[2],
                    "type": "scatter"
                }
                data_list2.append(insert_data2)
        self.create_db_table('bivariate_box_scatter', data_list2)

        update_dict1 = json.dumps({
                "category": "category",
                "name": "name",
                "group":"box_name",
                "condition": {'type': "box"}
            })

        update_dict2 = json.dumps({
                "category": "category",
                "data": ["x", "y"],
                "name": "sample",
                "condition": {'type': "scatter"}
            })
        self.update_db_record("bivariate_box", main_id, main_id=main_id, box_data=update_dict1, scatter_data=update_dict2, status="end")


