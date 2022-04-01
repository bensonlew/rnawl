# -*- coding: utf-8 -*-



import os
from bson.objectid import ObjectId
import pandas as pd
import glob
from api_base import ApiBase
import json
import numpy as np

class Pie(ApiBase):
    def __init__(self, bind_object):
        super(Pie, self).__init__(bind_object)

    def add_detail(self,main_id,  result_path):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        data_list = []

        data = pd.read_table(result_path, sep='\t',index_col=0)
        data.fillna(0,inplace=True)
        num_data = data.select_dtypes(include=[np.number])
        num_cols = num_data.columns
        for col in num_cols:
            tmp_sum = num_data[col].sum()
            num_data[col+'_per'] = num_data[col]/tmp_sum

        for index in data.index:
            for col in num_cols:
                insert_data = {
                    "pie_id" : main_id,
                    "row_name" : index,
                    "col_name" : col,
                    "value" : num_data[col][index],
                    "percent" : round(num_data[col+'_per'][index],5),
                    "type" : "pie"
                }
                data_list.append(insert_data)

        try:
            collection = self.db["pie_detail"]
            collection.insert_many(data_list)

            pie_data = json.dumps({"name":"row_name","data":"value","condition":{"type":"pie"}})
            col_name =  json.dumps({"data": num_cols.tolist()})
            self.db["pie"].update({"_id": main_id}, {"$set": {"pie_data":pie_data, "col_name":col_name}})

        except Exception as e:
            self.bind_object.logger.info("导入pie_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入pie_detail数据成功")

