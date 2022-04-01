# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
#  20210409

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import pandas as pd
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class O2plsda(ApiBase):
    def __init__(self, bind_object):
        super(O2plsda, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_o2plsda(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "O2plsda" + "_"
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='O2plsda',
                params=params,
                status="start")
            main_id = self.create_db_table('o2plsda', [main_info])
        else:
            main_id = ObjectId(main_id)

            try:
                self.update_db_record('o2plsda', main_id, )
            except Exception as e:
                self.bind_object.logger.error("导入o2plsda数据出错:%s" % e)
            else:
                self.bind_object.logger.info("导入o2plsda数据成功")
        return main_id

    def add_o2plsda_detail(self, main_id, file1, file2, file3, name):
        with open(file1, "r") as f:
            lines = f.readlines()
            data_list = []
            for line in lines[1:]:
                lin = line.strip().split("\t")
                insert_data = {
                    "o2pls_id": ObjectId(main_id),
                    "r2x": lin[0],
                    "r2y": lin[1],
                    "r2x_jonit": lin[2],
                    "r2y_jonit": lin[3],
                    "r2x_hat": lin[4],
                    "r2y_hat": lin[5],
                    "r2x_pred": lin[6],
                    "r2y_pred": lin[7],
                    "type": "table",
                }
            data_list.append(insert_data)
            self.create_db_table('o2plsda_stat', data_list)
        table_dict = {
            "column": [{"field": "r2x", "filter": "false", "sort": "false", "title": "R2X", "type": "float"},
            {"field": "r2y", "filter": "false", "sort": "false", "title": "R2Y", "type": "string"},
            {"field": "r2x_jonit","filter": "false","sort": "false","title": "R2X jonit","type": "string"},
            {"field": "r2y_jonit","filter": "false","sort": "false","title": "R2Y jonit","type": "string"},
            {"field": "r2x_hat","filter": "false","sort": "false","title": "R2X hat","type": "string"},
            {"field": "r2y_hat","filter": "false","sort": "false","title": "R2Y hat","type": "string"},
            {"field": "r2x_pred","filter": "false","sort": "false","title": "R2X pred","type": "string"},
            {"field": "r2y_pred","filter": "false","sort": "false","title": "R2Y pred","type": "string"}],
            "condition":{"type": "table"}}
        table1_info = json.dumps(table_dict, sort_keys=False, separators=(',', ':'))
        data2 = pd.read_table(file2, sep='\t', header=0)
        list1 = name.split(";")
        data_list1 = []
        print(list1)
        for i in list1:
            print(i)
            data1 = data2[data2["Omics"] == i]
            data1['sort'] = data1['pq[1]'].abs().astype(float)
            data1 = data1.sort_values(['sort'], ascending=False).drop('sort', axis=1)
            num =0
            for index, row in data1.iterrows():
                num +=1
                insert_data = {
                    "o2pls_id": ObjectId(main_id),
                    "name": row[0],
                    "x": row["pq[1]"],
                    "y": row["pq[2]"],
                    "sort": num,
                    "num": num,
                    "omics": row['Omics'],
                    'type': 'column',
                    }
                data_list1.append(insert_data)
        self.create_db_table('o2plsda_loading', data_list1)
        table2_dict = {"name":"name","data":["x", "y", "name"],"category":"omics", "condition":{"type":"scatter"}}
        table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
        table3_dict = {"name": "name", "data": "x", "category": "omics", "condition": {"type": "column"}}
        table3_info = json.dumps(table3_dict, sort_keys=False, separators=(',', ':'))
        data4 = pd.read_table(file3, sep='\t', header=0)
        data_list3 = []
        for index, row in data4.iterrows():
            insert_data = {
                "o2pls_id": ObjectId(main_id),
                "name": row["Sample"],
                "x": row["tu[1]"],
                "y": row["tu[2]"],
                "omics": row['Omics'],
                "category": row['Group'],
                'type':'scatter',
            }
            data_list3.append(insert_data)
        self.create_db_table('o2plsda_scores', data_list3)
        table4_dict = {"name": "name", "data": ["x", "y", "name"], "category": "omics", "category2": "omics", "condition": {"type": "scatter"}}
        table4_info = json.dumps(table4_dict, sort_keys=False, separators=(',', ':'))
        try:
            self.update_db_record('o2plsda', main_id, main_id=main_id, stat_table = table1_info, loading_data = table2_info, column_data = table3_info , scores_data=table4_info, status="end")
        except Exception as e:
            self.bind_object.logger.error("导入o2plsda_scores数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入o2plsda_scores数据成功")