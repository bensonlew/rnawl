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


class MirandaPic(ApiBase):
    def __init__(self, bind_object):
        super(MirandaPic, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_miranda_pic(self, pics=None, out_table=None, main_id=None, params=None, project_sn='miranda_pic', task_id='miranda_pic'):
        # add main table info
        if main_id is None:
            name = "miRanda_Pic" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='miRanda pic',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_miranda_pic', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        detail_data = {
            'column': [],
            'condition': {},
        }
        page_col = ['miRNA', 'Target', 'Total Score', 'Total Energy', 'Max Score', 'Max Energy', 'miRNA Length',
                    'Target Length', 'Position']
        data_col = ['mirna', 'target', 'total_score', 'total_energy', 'max_score', 'max_energy', 'mi_length',
                    'target_Length', 'position']
        length = len(data_col)
        for i in range(length):
            if i in [0, 1, length]:
                d_type = 'string'
            elif i in [6, 7]:
                d_type = 'int'
            else:
                d_type = 'float'
            col_detail = {
                "filter": "false",
                "field": data_col[i],
                "title": page_col[i],
                "type": d_type,
                "sort": "false",
            }
            detail_data['column'].append(col_detail)

        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'status': 'end',
            'main_id': main_id,
            'detail_data': json.dumps(detail_data, sort_keys=True, separators=(',', ':')),
            'hits': ['no'],
        }

        if len(out_table) > 0:
            self.add_miranda_detail(main_id=main_id, data=out_table[0])
            update_dict.update({'hits': ['yes']})
        self.update_db_record('sg_miranda_pic', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_miranda_detail(self, main_id, data):
        data_df = pd.read_table(data, header=0)
        data_df["miranda_id"] = main_id
        miranda_list = data_df.to_dict('records')
        try:
            self.col_insert_data('sg_miranda_pic_detail', miranda_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_miranda_pic_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入sg_miranda_pic_detail:%s成功" % (main_id,))