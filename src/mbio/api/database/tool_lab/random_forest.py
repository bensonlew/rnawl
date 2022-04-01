# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import types
import re
import datetime
import os
import json
from api_base import ApiBase

class RandomForest(ApiBase):
    def __init__(self, bind_object):
        super(RandomForest, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'
    
    def add_column_detail(self, main_id, input_file):
        with open(input_file,"r") as infile:
            line = infile.readline()
            # fd = line.rstrip().split('\t')
            # print(fd)
            # column_data_insert = []
            
            insert_dict = {
                    "condition" :{"type":"column"},
                    "data":"value",
                    "name":"name",
                }       
            # column_data_insert.append(insert_dict)
            update_dict ={
                "column_data":json.dumps(insert_dict)
            }
            self.update_db_record(
                "random_forest",{"_id":self.check_objectid(main_id)}, update_dict)
            n = 0
            while 1:
                insert_column_data = []
                line =infile.readline()
                if not line:
                    break
                n = n+1
                if n > 30 :
                    break
                fd = line.rstrip().split('\t')
                insert_column_dict = {
                    "random_forest_id": self.check_objectid(main_id),
                    "type":"column",
                    "name": fd[0],
                    "value":float(fd[1]),
                    "num": n
                }
                insert_column_data.append(insert_column_dict)
                self.col_insert_data("random_forest_column_detail",insert_column_data)
        insert_dict = {
                    "condition" :{"type":"line"},
                    "name":"name",
                }
            # line_data_insert.append(insert_dict)
        update_dict ={
                "line_data":json.dumps(insert_dict)
            }
        self.update_db_record(
                "random_forest",{"_id":self.check_objectid(main_id)}, update_dict)
        insert_dict = {
                    "condition" :{"type":"scatter"},
                    "data":["x","y"],
                    "category":"category",
                    "name":"name",
                }       
            # line_data_insert.append(insert_dict)
        update_dict ={
                "scatter_data":json.dumps(insert_dict)
            }
        self.update_db_record(
                "random_forest",{"_id":self.check_objectid(main_id)}, update_dict)



    def add_table_detail(self, main_id, input_file):
        with open(input_file,"r") as infile:
            line = infile.readline()
            # fd = line.rstrip().split('\t')
            # print(fd)
            fd = ["Feature","Importance"]
            table_data_insert =[]
            for i in fd:
                insert_dict = {
                    "field": i,
                    "title": i,
                    "sort":"false",
                    "filter":"false",
                    "type":"string",
                }
                table_data_insert.append(insert_dict)
            random_forest_table_dict = {
                "column" : table_data_insert,
                "condition":{'type':"table"}
                }
            
            update_dict = {
                "table_data": json.dumps(random_forest_table_dict)
            }
            self.update_db_record(
                "random_forest",{"_id":self.check_objectid(main_id)}, update_dict)
            
            while 1:
                insert_table = []
                line = infile.readline()
                if not line:
                    break
                fd1 = line.rstrip().split('\t')
                insert_table_dict = {
                    "type":"table",
                    "random_forest_id": self.check_objectid(main_id),
                }
                for i in range(len(fd)):
                    insert_table_dict[fd[i]] = fd1[i]
                    print("-----------")
                    print(fd[i])
                    print(fd1[i])
                    print("-"*10)
                insert_table.append(insert_table_dict)
                self.col_insert_data("random_forest_table_detail",insert_table)

    def add_line_detail(self, main_id, input_file):
        with open(input_file,"r") as infile:
            line = infile.readline()
            # fd = line.rstrip().split('\t')
            # print(fd)
            # line_data_insert = []
            while 1:
                insert_line_data = []
                insert_scatter_data = []
                line = infile.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                insert_line_dict = {
                    "random_forest_id": self.check_objectid(main_id),
                    "category":"group",
                    "type":"line",
                    "subtype" : "curveLinear",
                    "name":"line",
                    "x": float(fd[0]),
                    "y":float(fd[1])
                }
                insert_line_data.append(insert_line_dict)
                self.col_insert_data("random_forest_line_detail",insert_line_data)
                insert_scatter_dict = {
                    "random_forest_id": self.check_objectid(main_id),
                    "category":"group",
                    "type":"scatter",
                    "name":fd[0],
                    "x":float(fd[0]),
                    "y":float(fd[1])
                }
                insert_scatter_data.append(insert_scatter_dict)
                self.col_insert_data("random_forest_scatter_detail",insert_scatter_data)
                

    def add_tooltip(self,main_id, method):
        term_name = ""
        if method == "AUC":
            term_name = "AUC"
        elif method == "CV":
            term_name = "error rate"
        else:
            term_name = "none"
        group_data_insert = "Number of top important features"
        update_dict = {
            "x": group_data_insert
            }
        self.update_db_record("random_forest", {"_id": self.check_objectid(main_id)}, update_dict)
        group_data_insert =  term_name
        update_dict = {
            "y": group_data_insert
            }
        self.update_db_record("random_forest", {"_id": self.check_objectid(main_id)}, update_dict)

if __name__ == "__main__":
    a = RandomForest(None)
    main_id = "202104071635108545217538"
    important_path = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/random_forest/ten/RandomForest_importance_table.txt"
    a.add_column_detail(main_id,important_path)
    a.add_table_detail(main_id,important_path)
    a.add_line_detail(main_id,"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/random_forest/ten/RandomForest_10-fold_CV.txt")
    a.add_tooltip(main_id,"AUC")


