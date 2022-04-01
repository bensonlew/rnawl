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

class VcfStructure(ApiBase):
    def __init__(self,bind_object):
        super(VcfStructure, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_structrue(self, main_id,k, file):
        """
        导入structure堆叠图的数据（堆叠图）
        """
        with open(file,"r") as pop:
            while 1:
                insert_column = []
                line = pop.readline()
                if not line:
                    break
                fd = line.rstrip().split("\t")
                data = fd[1].split(' ')
                n = 0
                for i in data:
                    n += 1
                    insert_column_dict = {
                    "type": "column",
                    "vcf_structure_id":self.check_objectid(main_id),
                    "category":str(n),
                    "name": fd[0],
                    "group": int(k),
                    "value": float(i),
                }
                    insert_column.append(insert_column_dict)
                self.col_insert_data("vcf_structure_column",insert_column)
            
    
    def add_cverror_detail(self, main_id, cv_file):
        insert_dict1 = {
                "name":"name",
                "data":"value",
                "condition":{"type":"column"},
                # "group":"group",
            }
        update_dict1 = {
                "column_data": json.dumps(insert_dict1)
            }
        self.update_db_record(
                "vcf_structure",{"_id":self.check_objectid(main_id)},update_dict1
            )
        group_list = []
        # file_num = 0
        with open(cv_file,"r") as cvfile:
            while 1:
                insert_line_data = []
                insert_scatter_data = []
                line = cvfile.readline()
                if not line:
                    break
                # file_num += 1
                fd = line.rstrip().split('\t')
                group_list.append(float(fd[0]))
                insert_line_dict = {
                    "vcf_structure_id": self.check_objectid(main_id),
                    "category":"group",
                    "type":"line",
                    "subtype":"curveLinear",
                    "name":"line",
                    "x": float(fd[0]),
                    "y": float(fd[1]),
                }
                insert_line_data.append(insert_line_dict)
                self.col_insert_data("vcf_structure_cv_detial",insert_line_data)
                insert_scatter_dict = {
                    "vcf_structure_id": self.check_objectid(main_id),
                    "category":"group",
                    "type":"scatter",
                    "name":fd[0],
                    "x":float(fd[0]),
                    "y":float(fd[1])
                }
                insert_scatter_data.append(insert_scatter_dict)
                self.col_insert_data("vcf_structure_cv_detial",insert_scatter_data)
        insert_dict2 = {
                "data":group_list
            }
        update_dict2 = {
                "group": json.dumps(insert_dict2)
            }
        self.update_db_record(
                "vcf_structure",{"_id":self.check_objectid(main_id)},update_dict2
            )
        insert_dict = {
                    "condition" :{"type":"line"},
                    "name":"name",
                }
            # line_data_insert.append(insert_dict)
        update_dict ={
                "line_data":json.dumps(insert_dict)
            }
        self.update_db_record(
                "vcf_structure",{"_id":self.check_objectid(main_id)}, update_dict)
        insert_dict3 = {
                    "condition" :{"type":"scatter"},
                    "data":["x","y"],
                    "category":"category",
                    "name":"name",
                }       
            # line_data_insert.append(insert_dict)
        update_dict3 ={
                "scatter_data":json.dumps(insert_dict3)
            }
        self.update_db_record(
                "vcf_structure",{"_id":self.check_objectid(main_id)}, update_dict3)