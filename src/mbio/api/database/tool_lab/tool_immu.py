# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20201012
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

class ToolImmu(ApiBase):
    def __init__(self, bind_object):
        super(ToolImmu, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_immu(self, main_id, immu, method,params=None, project_sn='tool_lab', task_id='tool_lab',):
        # add main table info
        if main_id is None:
            name = "Immunedeconv" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='Immunedeconv',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_tool_immunedeconv', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(immu, header=0, sep='\t')
        immu_columns = df.columns
        columns_list = list()
        columns_list.append({'field': 'cell_type', 'filter': False, 'sort': False, 'title': 'cell_type', 'type': 'string'})
        for i in immu_columns[1:]:
            columns_list.append({'field': i, 'filter': False, 'sort': False, 'title': i, 'type': 'float'})
        data_columns = {'column': columns_list, 'condition': {}}
        columns_data = json.dumps(data_columns)

        # row_name = df['cell_type'].tolist()

        names = df.columns.tolist()[1:]
        # if method.lower() in ['xcell', 'cibersort']:
        #     cell_type = row_name[:-3]
        # else:
        #     cell_type = row_name
        df['immu_id'] = main_id
        detail = df.to_dict('r')
        self.create_db_table('sg_tool_immunedeconv_detail', detail)
        self.update_db_record('sg_tool_immunedeconv', main_id, status='end', main_id=main_id,
                              names=names, column_data_detail=columns_data)
        df = pd.read_table(immu, header=0, sep='\t')
        if method.lower() in ['xcell']:
            a = df.set_index('cell_type')
            a.drop(['immune score', 'stroma score', 'microenvironment score'], inplace=True)
            b = a.reset_index()
            b = b.T
        elif method.lower() in ['cibersort']:
            a = df.set_index('cell_type')
            a.drop(['P-value', 'Correlation', 'RMSE'], inplace=True)
            b = a.reset_index()
            b = b.T
        else:
            b = df.T
        c = list(b.iloc[0])
        d = b.drop(b.index[0])
        e = d.to_dict("r")
        sample = list(d.index)
        insert_data=[]
        for p,i in enumerate(e):
            for n, t in enumerate(c):
                single_detail=dict()
                single_detail["name"] = sample[p]
                single_detail["type"] = "column"
                single_detail['immu_id'] = main_id
                single_detail["category"] = t
                single_detail["value"] = i[n]
                insert_data.append(single_detail)
        try:
            self.create_db_table('sg_tool_immunedeconv_stacked_detail', insert_data)
            lista=list(df.columns[1:])
            dict_a = {"name": "name", "data": "value", "category": "category", "condition": {"type": "column"}}
            dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
            dict_c = {"data": lista}
            dict_d = json.dumps(dict_c, sort_keys=True, separators=(',', ':'))
            self.update_db_record('sg_tool_immunedeconv', main_id, status='end', main_id=main_id, column_data_stacked=dict_b,
                                  names=dict_d)
        except Exception as e:
            self.bind_object.logger.error("导入stacked_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入stacked_detail表格成功")
        if method.lower() == 'xcell':
            bar_df = df.set_index('cell_type').loc[['immune score', 'stroma score', 'microenvironment score']]
            b = bar_df.reset_index().T
            c = list(b.iloc[0])
            d = b.drop(b.index[0])
            e = d.to_dict("r")
            sample = list(d.index)
            insert_data = []
            for p, i in enumerate(e):
                for n, t in enumerate(c):
                    single_detail = dict()
                    single_detail["name"] = sample[p]
                    single_detail["type"] = "column"
                    single_detail['immu_id'] = main_id
                    single_detail["category"] = t
                    single_detail["value"] = i[n]
                    insert_data.append(single_detail)
            try:
                self.create_db_table('sg_tool_immunedeconv_bar_detail', insert_data)
                lista = list(df.columns[1:])
                dict_a = {"name": "name", "data": "value", "category": "category", "condition": {"type": "column"}}
                dict_b = json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
                dict_c = {"data": lista}
                dict_d = json.dumps(dict_c, sort_keys=True, separators=(',', ':'))
                self.update_db_record('sg_tool_immunedeconv', main_id, status='end', main_id=main_id, column_data_bar=dict_b,
                                      names_bar=dict_d)
            except Exception as e:
                self.bind_object.logger.error("导入stacked_detail表格信息出错:{}".format(e))
            else:
                self.bind_object.logger.info("导入stacked_detail表格成功")

