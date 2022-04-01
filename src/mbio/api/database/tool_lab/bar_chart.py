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


class BarChart(ApiBase):
    def __init__(self, bind_object):
        super(BarChart, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_bar_detail(self, data, main_id=None, params=None, data_file=None, group_file=None,
                       project_sn='bar_chart', task_id='bar_chart'):
        # add main table info
        if main_id is None:
            name = "Bar_Chart" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Bar Chart',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('bar_chart', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        self.add_each_bar(main_id=main_id, data=data)
        if group_file:
            bar_order = self.get_order(group_file=group_file)
        elif data_file:
            bar_order = self.get_order(data_file=data_file)

        query_dict = {
            '_id': main_id
        }

        update_dict = {
            "column_data": json.dumps({
                "name": "name",
                "data": "value",
                "category": "category",
                "condition": {'type': "column"}
            }),
            "ishape_data": json.dumps({
                "name": "name",
                "data": ["mean", "std_high", "std_low"],
                "group": "group",
                "condition": {'type': "ishape"}
            }),
            'bar_order': bar_order,
            'status': 'end',
            'main_id': main_id,
        }

        self.update_db_record('bar_chart', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_each_bar(self, main_id, data):
        data_df = pd.read_table(data, header=0)
        length = len(data_df.columns)
        ishape_list = list()
        data_df["bar_id"] = main_id
        data_df['category'] = data_df['name']
        if length == 3:
            ishape_df = data_df
            data_df = data_df.loc[:, ['name', 'mean', 'bar_id', 'category']]
            data_df.rename(columns={'mean': 'value'}, inplace=True)
            ishape_df["std_low"] = ishape_df['std']
            ishape_df.rename(columns={'std': 'std_high', 'category': 'group'}, inplace=True)
            ishape_df['type'] = 'ishape'
            ishape_list = ishape_df.to_dict('records')
        data_df['type'] = 'column'
        bar_list = data_df.to_dict('records')
        if ishape_list:
            bar_list.extend(ishape_list)
        try:
            self.col_insert_data('bar_chart_detail', bar_list)
        except Exception as e:
            self.bind_object.logger.info("导入bar_chart_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入bar_chart_detail:%s成功" % (main_id,))

    def get_order(self, data_file=None, group_file=None):
        if group_file:
            group_df = pd.read_table(group_file, header=None)
            bar_order = group_df.iloc[:, 1].drop_duplicates().values.tolist()
        elif data_file:
            data_df = pd.read_table(data_file, header=0)
            bar_order = data_df.iloc[:, 0].drop_duplicates().values.tolist()
        bar_order = ';'.join(bar_order)
        return bar_order
