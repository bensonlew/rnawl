# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

from bson.objectid import ObjectId
import datetime
import json
import os
import re
import pickle
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class LineChart(ApiBase):
    def __init__(self, bind_object):
        super(LineChart, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_line_chart(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "line_chart" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='line_chart',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('line_chart', [main_info])
        else:
            main_id = ObjectId(main_id)
            self.update_db_record('line_chart')
        return main_id

    def add_line_chart_detail(self, main_id, input_table, line_type=None, group_repetition=None):
        dfs = pd.read_csv(input_table, sep='\t')
        data_list = []
        data_list1 = []
        data_list2 = []
        group = []
        for i in dfs.index:
            insert_data = {"line_chart_id": main_id,  # 折线图详情表
                           "x": dfs.iloc[i, 0],
                           "type": "line"
                           }
            insert_data1 = {"line_chart_id": main_id,  # 散点详情表
                           "x": dfs.iloc[i, 0],
                           "type": "scatter"
                           }
            if group_repetition == "true":
                insert_data2 = {"line_chart_id": main_id,
                            "name": dfs.iloc[i, 0],
                            "type": "ishape"
                               }
            if line_type == "single":  # 单折线详情表
                insert_data.update({"name": "name",
                                    "category": "category",
                                    "y": float(dfs.iloc[i, 1])})
                data_list.append(insert_data)
                insert_data1.update({"name": "name",
                                     "category": "category",
                                     "y": float(dfs.iloc[i, 1])})
                data_list1.append(insert_data1)
                if group_repetition == "true":
                    insert_data2.update({"group": "all",
                                         "mean": float(dfs.iloc[i, 1]),
                                         "std_high": float(dfs.iloc[i, 2]),
                                         "std_low": float(dfs.iloc[i, 2])})
                    data_list2.append(insert_data2)
                group.append("group")
            else:  # 多折线详情表
                insert_data.update({"name": dfs.iloc[i, 1],
                                    "category": dfs.iloc[i, 1],
                                    "y": float(dfs.iloc[i, 2])})
                data_list.append(insert_data)
                insert_data1.update({"name": dfs.iloc[i, 1],
                                     "category": dfs.iloc[i, 1],
                                     "y": float(dfs.iloc[i, 2])})
                data_list1.append(insert_data1)
                if group_repetition == "true":
                    insert_data2.update({"group": dfs.iloc[i, 1],
                                         "mean": float(dfs.iloc[i, 2]),
                                         "std_high": float(dfs.iloc[i, 3]),
                                         "std_low": float(dfs.iloc[i, 3])})
                    data_list2.append(insert_data2)
                group.append(dfs.iloc[i, 1])
        try:
            self.create_db_table('line_chart_detail', data_list)
            self.create_db_table('line_chart_scatter_detail', data_list1)
            line_data_dict = {"name": "name", "category": "category" ,"condition": {"type": "line"}}
            line_data = json.dumps(line_data_dict, sort_keys=True, separators=(',', ':'))
            scatter_data_dict = {"name": "name", "category": "category", "data": ["x","y"], "condition": {"type": "scatter"}}
            scatter_data = json.dumps(scatter_data_dict, sort_keys=True, separators=(',', ':'))
            group_list = list(set(group))
            ishape_dict = {"name": "name", "data":["mean", "std_high", "std_low"], "condition": {'type': "ishape"}}
            ishape_data = json.dumps(ishape_dict, sort_keys=True, separators=(',', ':'))
            if group_repetition == "true":
                self.create_db_table('line_chart_ishape', data_list2)
                self.update_db_record('line_chart', main_id, status="end", main_id=main_id, line_data=line_data, scatter_data=scatter_data, ishape_data=ishape_data, ishape="yes", group=group_list)
            else:
                self.update_db_record('line_chart', main_id, status="end", main_id=main_id, line_data=line_data, scatter_data=scatter_data, ishape_data=ishape_data, ishape="no", group=group_list)
        except Exception as e:
            self.bind_object.logger.error("导入line_chart表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入line_chart表格成功")
