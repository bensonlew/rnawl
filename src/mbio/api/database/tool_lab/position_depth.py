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


class PositionDepth(ApiBase):
    def __init__(self, bind_object):
        super(PositionDepth, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_position_depth(self, result, main_id=None, params=None, project_sn='position_depth', task_id='position_depth'):
        # add main table info
        if main_id is None:
            name = "Position_Depth" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Position Depth',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_position_depth', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'status': 'end',
            'main_id': main_id,
        }
        if len(result) > 0:
            detail_data = self.add_depth_detail(main_id=main_id, data=result[0])
            update_dict.update({'detail_data': json.dumps(detail_data, sort_keys=True, separators=(',', ':'))})
        self.update_db_record('sg_position_depth', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_depth_detail(self, main_id, data):
        detail_data = {
            'column': [],
            'condition': {},
        }
        page_col = ['Chromosome', 'Position', 'Depth']
        data_df = pd.read_table(data, header=None, dtype=str)
        data_df.rename(columns={0: 'chr', 1: 'pos', 2: 'depth'}, inplace=True)
        length = len(data_df.columns)
        for i in range(length):
            d_type = 'string'
            col_detail = {
                "filter": "false",
                "field": data_df.columns[i],
                "title": page_col[i],
                "type": d_type,
                "sort": "false",
            }
            detail_data['column'].append(col_detail)
        data_df["depth_id"] = main_id
        depth_list = data_df.to_dict('records')
        try:
            self.col_insert_data('sg_position_depth_detail', depth_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_position_depth_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入sg_position_depth_detail:%s成功" % (main_id,))
        return detail_data