# -*- coding: utf-8 -*-

import os
from bson.objectid import ObjectId
import pandas as pd
import glob
from api_base import ApiBase
import json
import numpy as np
import copy

class Box(ApiBase):
    def __init__(self, bind_object):
        super(Box, self).__init__(bind_object)


    def add_detail(self,main_id,  result_path,group_path, method, log_change=None):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        data_list = []

        data = pd.read_table(result_path, sep='\t',index_col=0)
        data.fillna(0,inplace=True)
        if log_change=='log2':
            data = np.log2(data+1)
        elif log_change=='log10':
            data = np.log10(data+1)

        group = pd.read_table(group_path, sep='\t')


        group.rename(columns={group.columns[0]:'sample',group.columns[1]:'group'}, inplace=True)
        g_list = group['group'].tolist()
        g_sort = list()
        for i in g_list:
            if i not in g_sort:
                g_sort.append(i)

        s_list = group['sample'].tolist()
        s_sort = list()
        for i in s_list:
            if i not in s_sort:
                s_sort.append(i)

        g_map_s = dict(list(group.groupby(by='group')))

        #s_map_g = dict(list(group.groupby(by='sample')))
        s_map_g = dict(group.values.tolist())

        if method in ['mean', 'sum', 'median']:
            #for g in g_map_s.keys():
            for g in g_sort:
                if method=='mean':
                    data[g+'_merge'] = data[g_map_s[g]['sample'].tolist()].apply(lambda x: x.mean(),axis=1)
                elif method =='sum':
                    data[g+'_merge'] = data[g_map_s[g]['sample'].tolist()].apply(lambda x: x.sum(),axis=1)
                else:
                    data[g+'_merge'] = data[g_map_s[g]['sample'].tolist()].apply(lambda x: x.median(),axis=1)

            #for g in g_map_s.keys():
            for g in g_sort:
                res = data[g+'_merge'].describe()
                delq= res['75%'] - res['25%']

                ##异常点
                up_out_line = res['75%'] + delq*1.5
                down_out_line = res['25%'] - delq*1.5
                error_dot = data[g+'_merge'][data[g+'_merge'].map(lambda x: True if x > up_out_line or x < down_out_line else False)]
                normal_dot = data[g+'_merge'][data[g+'_merge'].map(lambda x: False if x > up_out_line or x < down_out_line else True)]
                nor_min = normal_dot.min()
                nor_max = normal_dot.max()

                insert_data = {
                    "box_id" : main_id,
                    "name" : g,
                    "q1" : res['25%'],
                    "q3" : res['75%'],
                    "median" : res['50%'],
                    "min" : nor_min,
                    "max" : nor_max,
                    "type" : "box"

                }

                data_list.append(insert_data)
                #tmp  = copy.deepcopy(insert_data)
                #tmp.update({'min': nor_min, 'max': nor_max, 'type': "box"})
                #data_list.append(tmp)

                ##异常点
                for index in error_dot.index:
                    scatter_data = {
                        "box_id" : main_id,
                        "name" : index,
                        "x" : g,
                        "y" : error_dot[index],
                        'type': "scatter",
                        'category': g

                    }
                    data_list.append(scatter_data)
        else:
            #for s in s_map_g.keys():
            for s in s_sort:
                res = data[s].describe()
                delq= res['75%'] - res['25%']

                ##异常点
                up_out_line = res['75%'] + delq*1.5
                down_out_line = res['25%'] - delq*1.5
                error_dot = data[s][data[s].map(lambda x: True if x > up_out_line or x < down_out_line else False)]
                normal_dot = data[s][data[s].map(lambda x: False if x > up_out_line or x < down_out_line else True)]
                nor_min = normal_dot.min()
                nor_max = normal_dot.max()

                insert_data = {
                    "box_id" : main_id,
                    "name" : s,
                    "category" : s_map_g[s],
                    "q1" : res['25%'],
                    "q3" : res['75%'],
                    "median" : res['50%'],
                    #"min" : res['min'],
                    #"max" : res['max'],
                    "min" : nor_min,
                    "max" : nor_max,
                    #"type" : "box1"
                    "type" : "box"
                }

                data_list.append(insert_data)
                # tmp  = copy.deepcopy(insert_data)
                # tmp.update({'min': nor_min, 'max': nor_max, 'type': "box"})
                # data_list.append(tmp)

                ##异常点
                for index in error_dot.index:
                    scatter_data = {
                        "box_id" : main_id,
                        "name" : index,
                        "x" : s,
                        "category" : s_map_g[s],
                        "y" : error_dot[index],
                        'type': "scatter"
                    }
                    data_list.append(scatter_data)

        try:
            collection = self.db["box_detail"]
            collection.insert_many(data_list)
            if method in ['mean', 'sum', 'median']:
                box_data = json.dumps({"name":"name","category":"name", "condition":{"type":"box"}})
            else:
                box_data = json.dumps({"name":"name","category":"category", "condition":{"type":"box"}})

            scatter_data = json.dumps({"name":"name", "data":["x", "y"],"category":"category","condition":{"type":"scatter"}})
            #box_data1 = json.dumps({"name":"name", "condition":{"type":"box1"}})
            if log_change=='log2':
                axis_str = 'log2(Value+1)'
            elif log_change=='log10':
                axis_str = 'log10(Value+1)'
            else:
                axis_str = 'Value'
            axis_data =  json.dumps({"data":[axis_str]})
            self.db["box"].update({"_id": main_id}, {"$set": {"box_data":box_data,"scatter_data":scatter_data, "main_id": main_id, "axis_data":axis_data}})
        except Exception as e:
            self.bind_object.logger.error("导入box_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入box_detail数据成功")


