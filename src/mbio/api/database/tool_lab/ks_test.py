# -*- coding: utf-8 -*-
# __author__ = 'zzg'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase
import os


class KsTest(ApiBase):
    def __init__(self, bind_object):
        super(KsTest, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_KsTest(self, file_path, project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
        if main_id is None:
            name = "KsTest"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='KsTest',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('sg_kstest', [main_info])
        else:
            main_id = ObjectId(main_id)

        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        
        table = [
            {"field": "sample_num", "type": "string", "sort": "false", "title": "sample_num"},
            {"field": "average", "type": "float", "sort": "false", "title": "average"},
            {"field": "std", "type": "float", "sort": "false", "title": "std"},
            {"field": "skewness", "type": "float", "sort": "false", "title": "skewness"},
            {"field": "kurtosis", "type": "float", "sort": "false", "title": "kurtosis"},
            {"field": "ks_d", "type": "float", "sort": "false", "title": "ks_D"},
            {"field": "ks_p", "type": "float", "sort": "false", "title": "ks_P"},
            {"field": "sw_d", "type": "float", "sort": "false", "title": "sw_D"},
            {"field": "sw_p", "type": "float", "sort": "false", "title": "sw_P"}
        ]
        table_data_update = {'column': table, 'condition': {}}
        column_data_insert = {
            "name": "name",
            "data": "value",
            "condition": {'type': "column"}
        }
        
        line_data_insert = {
            "name": "name",
            "x": "x",
            "y": "y",
            "condition": {'type': "line"}
        }

        scatter_data_insert = {
            "name": "name",
            "data": ["x","y"],
            "category": "category",
            "condition": {'type': "scatter"}
        }
        
        divide_line_data_insert = {
            "name": "name",
            "condition": {'type': "divide_line"}
        }

        table_info = json.dumps(table_data_update, sort_keys=True, separators=(',', ':'))
        column_data_info = json.dumps(column_data_insert, sort_keys=True, separators=(',', ':'))
        line_data_info = json.dumps(line_data_insert, sort_keys=True, separators=(',', ':'))
        scatter_data_info = json.dumps(scatter_data_insert, sort_keys=True, separators=(',', ':'))
        divide_line_data_info = json.dumps(divide_line_data_insert, sort_keys=True, separators=(',', ':'))
        
        out1_file = os.path.join(file_path, 'out1.xlsx')
        histogram_X_file = os.path.join(file_path, 'histogram_X_axis.txt')
        histogram_Y_file = os.path.join(file_path, 'histogram_Y_axis.txt')
        histogram_line_X = os.path.join(file_path, 'histogram_line_X.txt')
        histogram_line_Y = os.path.join(file_path, 'histogram_line_Y.txt')
        QQ_line_file = os.path.join(file_path, 'QQ_line.txt')
        QQ_X_file = os.path.join(file_path, 'QQ_X_axis.txt')
        QQ_Y_file = os.path.join(file_path, 'QQ_Y_axis.txt')
        insert_data1 = [];insert_data2 = [];insert_data3 = [];insert_data4 = [];insert_data5 = [];insert_data6 = []
        data1 = [];data2 = [];data3 = [];data4 = [];data7 = [];data8 = []

        #table_detil
        with open(out1_file, 'r') as r5:
            lines5 = r5.readlines()
            lines5_column2 = lines5[1].strip().split("\t")
            insert_data3.append({
                "sample_num":lines5_column2[0],
                "average":lines5_column2[1],
                "std":lines5_column2[2],
                "skewness":lines5_column2[3],
                "kurtosis":lines5_column2[4],
                "ks_d":lines5_column2[5],
                "ks_p":lines5_column2[6],
                "sw_d":lines5_column2[7],
                "sw_p":lines5_column2[8],
                "kstest_id": main_id
            })
        self.create_db_table("sg_kstest_table", insert_data3)
        
        #scatter_detil
        with open(QQ_X_file, 'r') as r1:
            with open(QQ_Y_file, 'r') as r2:
                lines1 = r1.readlines()
                lines2 = r2.readlines()
                for line1 in lines1:
                    data1.append(line1.strip())
                for line2 in lines2:
                    data2.append(line2.strip())
                for i in range(len(data1)):
                    insert_data1.append({
                        "category": "All",
                        "kstest_id": main_id,
                        "x":str(round(float((data1[i])))),
                        "y":float(data2[i]),
                        "type":"scatter"
                    })
        self.create_db_table("sg_kstest_scatter", insert_data1)
            
        #line_detil
        with open(histogram_line_X, 'r') as r7:
            with open(histogram_line_Y, 'r') as r8:
                lines7 = r7.readlines()
                lines8 = r8.readlines()
                for line7 in lines7:
                    data7.append(line7.strip())
                for line8 in lines8:
                    data8.append(line8.strip())
                for i in range(len(data8)):
                    insert_data6.append({
                        "kstest_id": main_id,
                        "x":round(float(data7[i])),
                        "y":float(data8[i]),
                        "type":"line"
                    })
        self.create_db_table("sg_kstest_line", insert_data6)
        
        #column_detil
        with open(histogram_X_file, 'r') as r3:
            with open(histogram_Y_file, 'r') as r4:
                lines3 = r3.readlines()
                lines4 = r4.readlines()
                for line3 in lines3:
                    data3.append(line3.strip())
                for line4 in lines4:
                    data4.append(line4.strip())                
                for j in range(len(data4)):
                    insert_data2.append({
                        "name": round(float(data3[j])),
                        "value": float(data4[j]),
                        "kstest_id": main_id,
                        "type": "column"
                    })
        self.create_db_table("sg_kstest_column", insert_data2)
        
        #divide_line_detil
        with open(QQ_line_file, 'r') as r6:
            line6 = r6.readlines()
            spot1 = line6[0].strip().split("\t")
            spot2 = line6[1].strip().split("\t")
            insert_data5.append({
                "name":"divide_line",
                "x1": spot1[0][3:],
                "y1": spot1[1][3:],
                "x2": spot2[0][3:],
                "y2": spot2[1][3:],
                "kstest_id": main_id,
                "type": "divide_line"
                    })
        self.create_db_table("sg_kstest_divide_line", insert_data5)

        self.update_db_record('sg_kstest', main_id, status="end", main_id=main_id, table_data=table_info, column_data=column_data_info, line_data=line_data_info, scatter_data=scatter_data_info, divide_line_data= divide_line_data_info)
        return main_id