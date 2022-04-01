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


class Bsa(ApiBase):
    def __init__(self, bind_object):
        super(Bsa, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_bsa_main(self, result=None, png=None, main_id=None, params=None, project_sn='bsa', task_id='bsa'):
        # add main table info
        if main_id is None:
            name = "BSA" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='BSA',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_bsa', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        detail_data = {
            'column': [],
            'condition': {},
        }
        page_col = ['Chromosome', 'Start', 'End', 'SNP Number']
        data_col = ['chr', 'start', 'end', 'num']
        length = len(data_col)
        for i in range(length):
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

        if png:
            pic = png
        else:
            pic = ''

        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'status': 'end',
            'main_id': main_id,
            'detail_data': json.dumps(detail_data, sort_keys=True, separators=(',', ':')),
            'pic_data': pic
        }
        if result and os.path.getsize(result) > 0:
            self.add_bsa_detail(main_id=main_id, data=result)
        self.update_db_record('sg_bsa', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_bsa_detail(self, main_id, data):
        data_df = pd.read_table(data, header=None, dtype={0: 'str'})
        data_df.rename(columns={0: 'chr', 1: 'start', 2: 'end', 3: 'num'}, inplace=True)
        data_df["bsa_id"] = main_id
        bsa_list = data_df.to_dict('records')
        try:
            self.col_insert_data('sg_bsa_detail', bsa_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_bsa_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入sg_bsa_detail:%s成功" % (main_id,))