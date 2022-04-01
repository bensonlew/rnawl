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

class Anosim(ApiBase):
    def __init__(self, bind_object):
        super(Anosim, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_box_detail(self, main_id, input_file,method):
        """
        导入箱线图数据
        """
        with open(input_file,"r") as infile:
            infile.readline()
            # fd = line.rstrip().split('\t')
            # print(fd)
            # box_data_insert =[]
            insert_dict = {
                    "condition" :{"type":"box"},
                    "name":"name",
                }       
            # box_data_insert.append(insert_dict)
            update_dict = {
                "box_data": json.dumps(insert_dict)
            }
            self.update_db_record(
                "anosim",{"_id":self.check_objectid(main_id)}, update_dict)
            
            while 1:
                insert_box = []
                line = infile.readline()
                if not line:
                    break
                fd1 = line.rstrip().split('\t')
                insert_table_dict = {
                    "type":"box",
                    "anosim_id": self.check_objectid(main_id),
                    "category":fd1[0],
                    "name":fd1[0],
                    "min":float(fd1[1]),
                    "q1":float(fd1[2]),
                    "median":float(fd1[4]),
                    "q3":float(fd1[5]),
                    "max":float(fd1[6]),
                }
                
                insert_box.append(insert_table_dict)
                self.col_insert_data("anosim_box_detail",insert_box)
        update_dict = {
            "method": method
        }
        self.update_db_record(
            "anosim",{"_id":self.check_objectid(main_id)}, update_dict)
    
    def add_table_detail(self,main_id,input_file):
        """
        导入表格数据
        """
        with open(input_file,"r") as infile:
            line = infile.readline()
            fd = line.rstrip().split('\t')
            print(fd)
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
            anosim_table_dict = {
                "column" : table_data_insert,
                "condition":{'type':"table"}
                }
            
            update_dict = {
                "table_data": json.dumps(anosim_table_dict)
            }
            self.update_db_record(
                "anosim",{"_id":self.check_objectid(main_id)}, update_dict)
            
            while 1:
                insert_table = []
                line = infile.readline()
                if not line:
                    break
                fd1 = line.rstrip().split('\t')
                insert_table_dict = {
                    "type":"table",
                    "anosim_id": self.check_objectid(main_id),
                }
                for i in range(len(fd)):
                    
                    insert_table_dict[fd[i].split('.')[0]] = fd1[i]
                    print("-----------")
                    print(fd[i])
                    print(fd1[i])
                    print("-"*10)
                insert_table.append(insert_table_dict)
                self.col_insert_data("anosim_table_detail",insert_table)


if __name__ == "__main__":
    a = Anosim(None)
    a.add_box_detail("202104071029108545217540","/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/anosim/box_data.xls","bray_curtis")
    a.add_table_detail("202104071029108545217540","/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/anosim/anosim_results.txt")