# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

from bson.objectid import ObjectId
import datetime
import json
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class ClusterHistogram(ApiBase):
    def __init__(self, bind_object):
        super(ClusterHistogram, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_cluster_histogram(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "cluster_histogram" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='cluster_histogram',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('cluster_histogram', [main_info])
        else:
            main_id = ObjectId(main_id)
            self.update_db_record('cluster_histogram')
        return main_id

    def add_cluster_histogram_detail(self, main_id, input_table):  # 柱状图数据详情表
        dfs = pd.read_csv(input_table, sep='\t')
        data_list = []
        for i in dfs.index:
            insert_data = {
                "cluster_histogram_id": main_id,
                "category": dfs.iloc[i, 1],
                "name": dfs.iloc[i, 0],
                "value": float(dfs.iloc[i, 2]),
                "type": "column"
            }
            data_list.append(insert_data)
        try:
            self.create_db_table('cluster_histogram_detail', data_list)
            column_data_dict = {"name": "name", "data": "value", "category": "category",
                                "condition": {"type": "column"}}
            column_data = json.dumps(column_data_dict, sort_keys=True, separators=(',', ':'))
            ishape_dict = {"name": "name", "data":["mean", "std_high", "std_low"], "condition": {'type': "ishape"}}
            ishape_data = json.dumps(ishape_dict, sort_keys=True, separators=(',', ':'))
            self.update_db_record('cluster_histogram', main_id, status="end", main_id=main_id,
                                  colunmn_data=column_data, ishape_data=ishape_data)
                                  
        except Exception as e:
            self.bind_object.logger.error("导入cluster_histogram_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入cluster_histogram_detail表格成功")

    def add_ishape_detail(self, main_id, result_file, error_bar=None):  # 工字数据详情表
        df_result = pd.read_csv(result_file, sep='\t')
        ishape_data_list = []
        if error_bar is "std":
            for i in df_result.index:
                insert_ishape_data = {
                    "cluster_histogram_id": main_id,
                    "group": df_result.iloc[i, 1],
                    "name": df_result.iloc[i, 0],
                    "mean": float(df_result.iloc[i, 2]),
                    "std_high": float(df_result.iloc[i, 3]),
                    "std_low": float(df_result.iloc[i, 3]),
                    "type": "ishape"
                }
                ishape_data_list.append(insert_ishape_data)
            try:
                self.create_db_table('cluster_histogram_ishape', ishape_data_list)
                ishape_data_dict = {"name": "name", "data": ["mean", "std_high", "std_low"],
                                    "condition": {'type': "ishape"}}
                ishape_data = json.dumps(ishape_data_dict, sort_keys=True, separators=(',', ':'))
                self.update_db_record('cluster_histogram', main_id, status="end", main_id=main_id,
                                      ishape_data=ishape_data)
            except Exception as e:
                self.bind_object.logger.error("导入cluster_histogram_ishape表格信息出错:{}".format(e))
            else:
                self.bind_object.logger.info("导入cluster_histogram_ishape表格成功")
        else:
            for i in df_result.index:
                insert_ishape_data = {
                    "cluster_histogram_id": main_id,
                    "group": df_result.ix[i, 1],
                    "name": df_result.ix[i, 0],
                    "mean": float(df_result.ix[i, 2]),
                    "sem_high": float(df_result.iloc[i, 3]),
                    "sem_low": float(df_result.iloc[i, 3]),
                    "type": "ishape"
                }
                ishape_data_list.append(insert_ishape_data)
            try:
                self.create_db_table('cluster_histogram_ishape', ishape_data_list)
                ishape_data_dict = {"name": "name", "data": ["mean", "sem_high", "sem_low"],
                                    "condition": {'type': "ishape"}}
                ishape_data = json.dumps(ishape_data_dict, sort_keys=True, separators=(',', ':'))
                self.update_db_record('cluster_histogram', main_id, status="end", main_id=main_id,
                                      ishape_data=ishape_data)
            except Exception as e:
                self.bind_object.logger.error("导入cluster_histogram_ishape表格信息出错:{}".format(e))
            else:
                self.bind_object.logger.info("导入cluster_histogram_ishape表格成功")
