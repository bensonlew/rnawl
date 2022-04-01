# !/usr/bin/python
# -*- coding: utf-8 -*-
from bson.objectid import ObjectId
from api_base import ApiBase
import json
import pandas as pd
import datetime
import glob,re
from mbio.api.database.ref_rna_v2.api_base import ApiBase


class AllCircle(ApiBase):
    def __init__(self, bind_object):
        super(AllCircle, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_circle(self, main_id, circ_file):
        df = ''
        with open(circ_file) as f:
            if re.search("\t",f.readlines()[0].strip()):
                df = pd.read_table(circ_file, index_col=0)
            else:
                df = pd.read_excel(circ_file, index_col=0)
        data1 = [{'name': i, 'value': 0, 'relation': [], 'sub_type': 'circle_right', 'type': 'cluster_chord', 'circle_id': ObjectId(main_id)} for i in df.columns.tolist()]
        data2 = [{'name': i, 'value': 0, 'relation': df.loc[i].tolist(), 'sub_type': 'circle_left', 'type': 'cluster_chord', 'circle_id': ObjectId(main_id)} for i in df.index.tolist()]
        self.create_db_table('sg_all_circle_detail', data1)
        self.create_db_table('sg_all_circle_detail', data2)
        dict_a = {"name": "name", "value": "value", "relation": "relation", "condition": {"type": "cluster_chord"}}
        dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_all_circle', main_id, cluster_chord_data=dict_b, status='end')

