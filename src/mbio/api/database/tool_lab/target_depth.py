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


class TargetDepth(ApiBase):
    def __init__(self, bind_object):
        super(TargetDepth, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_target_depth(self, out_table, seq_len, main_id=None, params=None, project_sn='target_depth', task_id='target_depth'):
        # add main table info
        if main_id is None:
            name = "Target_Depth" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Target Depth',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_target_depth', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        detail_data = self.add_target_detail(main_id=main_id, data=out_table, seq_len=seq_len)
        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'status': 'end',
            'main_id': main_id,
            'detail_data': json.dumps(detail_data, sort_keys=True, separators=(',', ':')),
        }
        self.update_db_record('sg_target_depth', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_target_detail(self, main_id, data, seq_len):
        depth_list = list()
        for each in data:
            sample = os.path.basename(each).split('_map.xls')[0]
            data_df = pd.read_table(each, header=0)
            total_depth = sum(data_df.iloc[:, 2])
            depth_data = {
                'sample': sample,
                'ave_depth': int(round(total_depth/float(seq_len))),
                'depth_id': main_id,
            }
            depth_list.append(depth_data)

        detail_data = {
            'column': [],
            'condition': {},
        }
        page_col = ['Sample Name', 'Ave Depth']
        data_col = ['sample', 'ave_depth']
        for i in range(len(data_col)):
            if i == 0:
                d_type = 'string'
            else:
                d_type = 'int'
            col_detail = {
                "filter": "false",
                "field": data_col[i],
                "title": page_col[i],
                "type": d_type,
                "sort": "false",
            }
            detail_data['column'].append(col_detail)
        try:
            self.col_insert_data('sg_target_depth_detail', depth_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_target_depth_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入sg_target_depth_detail:%s成功" % (main_id,))
        return detail_data
