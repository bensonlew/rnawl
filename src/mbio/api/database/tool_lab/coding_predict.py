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


class CodingPredict(ApiBase):
    def __init__(self, bind_object):
        super(CodingPredict, self).__init__(bind_object)
        self._project_type = 'tool_lab'


    def add_coding_stat_detail  (self, predict_results,main_id=None):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(predict_results):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(predict_results))

        if main_id is None:
            name = "only_for_test"
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                # project_sn=project_sn,
                task_id="coding_predict",
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
            main_id = self.create_db_table('coding_predict', [main_info])

        df = pd.read_table(predict_results,header=None)
        df = df.fillna(0)
        df.columns = ["name", "seqs"]
        detail_dict_list = df.to_dict('records')
        for i in detail_dict_list:
            if i["seqs"]  == 0:
                i["seqs"] = []
        self.create_db_table('coding_predict_detail', detail_dict_list, tag_dict={'coding_predict_id': main_id, 'type': 'venn'})
        venn_data = dict(names='seqs', category='name', condition={'type': 'venn'})
        venn_data = json.dumps(venn_data, sort_keys=True, separators=(',', ':'))
        self.update_db_record('coding_predict', main_id, status="end", main_id=main_id, venn_data=venn_data)

