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


class FastqDup(ApiBase):
    def __init__(self, bind_object):
        super(FastqDup, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_fastq_dup(self, data, main_id=None, params=None, project_sn='fastq_dup', task_id='fastq_dup'):
        # add main table info
        if main_id is None:
            name = "Fastq_Dup" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Fastq dup',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_fastq_dup', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        if os.path.exists(data):
            detail_data = self.add_fastq_dup_detail(main_id=main_id, data=data)
            query_dict = {
                '_id': main_id
            }

            update_dict = {
                'detail_data': json.dumps(detail_data, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'main_id': main_id,
            }
            self.update_db_record('sg_fastq_dup', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_fastq_dup_detail(self, main_id, data):
        detail_data = {
            'column': [],
            'condition': {},
        }
        page_col = ['Sample', 'Dup_R1(%)', 'Dup_R2(%)', 'Dup_Pair(%)']
        data_df = pd.read_table(data, header=0)
        length = len(data_df.columns)
        if length == 4:
            result_df = data_df.rename(columns={'Sample': 'sample', 'Dup_R1(%)': 'dup_r1', 'Dup_R2(%)': 'dup_r2',
                                                'Dup_Pair(%)': 'dup_pair'})
        else:
            result_df = data_df.rename(columns={'Sample': 'sample', 'Dup_R1(%)': 'dup_r1'})

        for i in range(len(result_df.columns)):
            if i == 0:
                type = 'string'
            else:
                type = 'float'
            col_detail = {
                "filter": "false",
                "field": result_df.columns[i],
                "title": page_col[i],
                "type": type,
                "sort": "false",
            }
            detail_data['column'].append(col_detail)

        result_df["dup_id"] = main_id
        dup_list = result_df.to_dict('records')
        try:
            self.col_insert_data('sg_fastq_dup_detail', dup_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_fastq_dup_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入sg_fastq_dup_detail:%s成功" % (main_id,))
        return detail_data
