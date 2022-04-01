# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

from bson.objectid import ObjectId
import datetime
import json
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class GroupScatter(ApiBase):
    def __init__(self, bind_object):
        super(GroupScatter, self).__init__(bind_object)
        self._project_type ='tool_lab'

    def add_group_scatter(self,project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
        if main_id is None:
            name = "group_scatter" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='group_scatter',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('group_scatter', [main_info])
        else:
            main_id = ObjectId(main_id)
            self.update_db_record('stacked_column_absolute')
        return main_id
    
    def add_group_scatter_detail(self, data_file, group_file, main_id):
        df_data = pd.read_csv(data_file,sep='\t')
        df_group = pd.read_csv(group_file,sep='\t')
        df_data.columns = ['sample_name','value']
        df_group.columns = ['sample_name','group']
        dfs = pd.merge(df_group, df_data,on='sample_name')
        data_list = []
        for j in dfs.index:
            insert_data_scatter = {
                "group_scatter_id": main_id,
                "category":dfs.iloc[j,1],
                "name": dfs.iloc[j,0],
                "x": dfs.iloc[j,1],
                "y": float(dfs.iloc[j,2]),
                'type': "scatter"
                }
            data_list.append(insert_data_scatter)
        try:
            collection = self.db['group_scatter_detail']
            collection.insert_many(data_list)
            scatter_data = json.dumps({ "category":"category", "name":"name", "data":["x", "y"],"condition":{"type":"scatter"}}, sort_keys=True, separators=(',',':'))
            self.update_db_record('group_scatter',main_id , status = "end",main_id=main_id, scatter_data=scatter_data)
            self.bind_object.logger.info("更新主表完成")
        except Exception as e:
            self.bind_object.logger.error("导入scatter_group_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入scatter_group_detail表格成功")

    def add_group_scatter_box_detail(self, result_file, main_id):
        df_result = pd.read_csv(result_file, sep='\t')
        data_list = []
        for i in df_result.index:
            insert_data = {
                "group_scatter_id": main_id,
                "name": df_result.iloc[i,0],
                "category": df_result.iloc[i,0],
                "q1": float(df_result.iloc[i,2]),
                "q3": float(df_result.iloc[i,4]),
                "median": float(df_result.iloc[i,3]),
                "min": float(df_result.iloc[i,1]),
                "max": float(df_result.iloc[i,5]),
                "type": "box"
                }
            data_list.append(insert_data)
        try:
            collection = self.db['group_scatter_box']
            collection.insert_many(data_list)
            box_data = json.dumps({"name": "name", "condition": {"type": "box"}},sort_keys=True,separators=(',',':'))
            self.update_db_record('group_scatter', main_id, status="end", main_id=main_id, box_data=box_data)
            self.bind_object.logger.info("更新主表完成")
        except Exception as e:
            self.bind_object.logger.error("导入scatter_group_box表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入scatter_group_box表格成功")

