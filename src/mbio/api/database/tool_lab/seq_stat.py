# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from bson.objectid import ObjectId
import types
import os


class SeqStat(ApiBase):
    def __init__(self, bind_object):
        super(SeqStat, self).__init__(bind_object)
        self._project_type = 'tool_lab'


    def add_seq_stat_detail(self, seq_stat_result,seq_length_result,main_id=None):
        # if not isinstance(main_id, ObjectId):
        #     if isinstance(main_id, types.StringTypes):
        #         main_table_id = ObjectId(main_id)
        #     else:
        #         self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(seq_stat_result):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(seq_stat_result))

        if main_id is None:
            name = "only_for_test"
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                # project_sn=project_sn,
                task_id="seq_stat",
                name=name,
                # exp_level=exp_level,
                # exp_type=exp_type.upper(),
                # method=quant_method,
                desc='stacked main table',
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                # params=params,
                version="v1",
                # sample_order=sample_order + group_order,
                status="start"
            )
            main_id = self.create_db_table('seq_stat', [main_info])

        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_table_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！')

        df = pd.read_table(seq_stat_result,header=None)
        df.columns=["range","num"]
        # df_row = list(df[df_columns[0]])
        # df_col = df_columns[1:]
        df['seq_stat_id'] = main_id
        df['type'] = "columns"
        df["range"] = df["range"].apply(lambda x: "_".join([str(int(float(i))) for i in x.split("-")]))
        column_dict={"name":"range","data":"num","condition":{"type":"column"}}
        columns_info = json.dumps(column_dict, sort_keys=True, separators=(',', ':'))
        table_dict = {"column":[{"field":"range","title":"SeqLength(bp)","sort":"false","type":"string","filter":"false"},
                                {"field": "num", "title": "SeqNum", "sort": "false", "type": "int","filter": "false"}],
                      "condition":{"type":"table"}}
        # table_dict={"range":"SeqLength(bp)","num":"SeqNum"}
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
        df2 = pd.read_table(seq_length_result, header=None)
        df2.columns = ["seq_name", "seq_length"]
        df2['seq_stat_id'] = main_id
        table_dict2 = {"column": [
            {"field": "seq_name", "title": "SeqName", "sort": "false", "type": "string", "filter": "false"},
            {"field": "seq_length", "title": "SeqLength(bp)", "sort": "true", "type": "int", "filter": "false"}],
                      "condition": {"type": "table"}}
        # table_dict2 = {"seq_name": "SeqName", "seq_length": "SeqLength(bp)"}
        table_info2 = json.dumps(table_dict2, sort_keys=True, separators=(',', ':'))
        self.create_db_table('seq_stat_detail', df2.to_dict('r'))
        self.create_db_table('seq_stat_detail', df.to_dict('r'))
        self.update_db_record('seq_stat', main_id, insert_dict={'main_id': main_id, 'status': 'end',
                                                                'table_data':table_info,'column_data':columns_info,'table_data_detail': table_info2})

