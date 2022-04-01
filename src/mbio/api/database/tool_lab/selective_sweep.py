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


class SelectiveSweep(ApiBase):
    def __init__(self, bind_object):
        super(SelectiveSweep, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_selective_sweep(self, result, png, pdf, main_id=None, params=None, project_sn='selective_sweep', task_id='selective_sweep'):
        # add main table info
        if main_id is None:
            name = "Selective_Sweep" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Selective Sweep',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_selective_sweep', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        detail_data = {
            'column': [],
            'condition': {},
        }
        page_col = ['Chromosome', 'Position', 'Pi1', 'Pi2', 'TajimaD1', 'TajimaD2', 'Fst']
        data_col = ['chr', 'pos', 'pi1', 'pi2', 'tajima1', 'tajima2', 'fst']
        length = len(data_col)
        for i in range(length):
            if i == 0:
                d_type = 'string'
            elif i == 1:
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
            'pic_data': png
        }

        if os.path.exists(result):
            self.add_sweep_detail(main_id=main_id, data=result)
        self.update_db_record('sg_selective_sweep', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_sweep_detail(self, main_id, data):
        data_df = pd.read_table(data, header=0, dtype={'chr': 'str', 'pos': 'int64'})
        data_df["sweep_id"] = main_id
        sweep_list = data_df.to_dict('records')
        try:
            self.col_insert_data('sg_selective_sweep_detail', sweep_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_selective_sweep_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入sg_selective_sweep_detail:%s成功" % (main_id,))