# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20210202
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
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class ProteinPicture(ApiBase):
    def __init__(self, bind_object):
        super(ProteinPicture, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_protein_picture(self, main_id, df_table, params=None, project_sn='tool_lab', task_id='tool_lab'):
        if main_id is None:
            name = "Protein_from_picture" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Protein from picture',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_protein_from_picture', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(df_table, header=0, sep='\t')
        df.rename(columns={'Total Spectrum': 'total_spectrum', 'Identified Spectrum': 'identified_spectrum',
                           'Peptide number': 'peptide_number', 'Protein number': 'protein_number',
                           'Protein group number': 'protein_group_number'}, inplace=True)
        columns_list = list()
        columns_list.append({'field': 'total_spectrum', 'filter': False, 'sort': False, 'title': 'Total Spectrum', 'type': 'int'})
        columns_list.append({'field': 'identified_spectrum', 'filter': False, 'sort': False, 'title': 'Identified Spectrum', 'type': 'int'})
        columns_list.append({'field': 'peptide_number', 'filter': False, 'sort': False, 'title': 'Peptide number', 'type': 'int'})
        columns_list.append({'field': 'protein_number', 'filter': False, 'sort': False, 'title': 'Total Spectrum', 'type': 'int'})
        columns_list.append({'field': 'protein_group_number', 'filter': False, 'sort': False, 'title': 'Identified Spectrum', 'type': 'int'})
        data_columns = {'column': columns_list, 'condition': {}}
        columns_data = json.dumps(data_columns)
        df['picture_id'] = main_id
        detail = df.to_dict('r')
        self.create_db_table('sg_protein_from_picture_detail', detail)
        self.update_db_record('sg_protein_from_picture', main_id, main_id=main_id, column_data_detail=columns_data)
        info = pd.read_table(df_table, header=None, index_col=None, sep='\t')
        info_new = info.T
        info_new.rename(columns={0:'name',1:'value'}, inplace=True)
        info_new['value'] = info_new['value'].astype(int)
        info_new["type"] = "column"
        info_new['picture_id'] = main_id
        info_list_detail = info_new.to_dict('r')
        self.create_db_table('sg_protein_from_picture_column_detail', info_list_detail)
        dict_a = {"name": "name", "data": "value", "condition": {"type": "column"}}
        dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_protein_from_picture', main_id, status='end', main_id=main_id, column_data_column=dict_b)



