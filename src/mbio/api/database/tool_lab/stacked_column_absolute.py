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


class StackedColumnAbsolute(ApiBase):
    def __init__(self, bind_object):
        super(StackedColumnAbsolute, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_stacked_column_absolute(self,project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
        if main_id is None:
            name = "stacked_column_absolute" + "_"
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name = name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='stacked_column_absolute',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('stacked_column_absolute', [main_info])
        else:
            main_id = ObjectId(main_id)
            self.update_db_record('stacked_column_absolute')
        return main_id

    def add_stacked_column_absolute_detail(self,  result_file, main_id):
        df_result = pd.read_csv(result_file,sep='\t')
        column_list = df_result.columns.tolist()
        data_list = []
        for i in df_result.index:
            for j in column_list[1:]:
                insert_data = {
                    "stacked_column_absolute_id": main_id,
                    "name": j,
                    "category": df_result.ix[i,0],
                    "value":float(df_result.ix[i,j]),
                    "type": "column"
                }
                data_list.append(insert_data)
        try:
            #self.create_db_table('stacked_column_absolute_detail', data_list)
            collection = self.db['stacked_column_absolute_detail']
            collection.insert_many(data_list)
            column_data_dict = {"name":"name","data":"value","category":"category","condition":{"type":"column"}}
            column_data = json.dumps(column_data_dict,sort_keys=True,separators=(',',':'))
            self.update_db_record('stacked_column_absolute', main_id, status="end", main_id=main_id, colunmn_data=column_data)
            self.bind_object.logger.info("更新主表完成")
        except Exception as e:
            self.bind_object.logger.error("导入stacked_column_absolute_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入stacked_column_absolute_detail表格成功")

