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

class Vpa(ApiBase):
    def __init__(self, bind_object):
        super(Vpa, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_detail(self, main_id, input_file):
        """
        导入VennVPA图数据
        """
        with open(input_file, "r") as infile:
            infile.readline()
            self.Residuals = 0.0
            while 1:
                insert_venn = []
                line = infile.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                if fd[0] == "Residuals":
                    if fd[1] != "NA":
                        self.Residuals = round(float(fd[1]),6)
                    else:
                        self.Residuals = fd[1]
                    break
                insert_data = fd[0].split('&')
                insert_table_dict = {
                    "type" : "venn",
                    "vpa_id":self.check_objectid(main_id),
                    "name" : fd[0],
                    "data" : insert_data,
                    "value": round(float(fd[1]),6) if fd[1] != "NA" else fd[1],
                    "category": fd[0]
                }
                insert_venn.append(insert_table_dict)
                self.col_insert_data("vpa_analysis_detail",insert_venn)
        insert_dict = {
                "names":"data",
                "value":"value",
                'category':"name",
                "condition": {"type":"venn"}
        }
        update_dict = {
                "venn_data": json.dumps(insert_dict)}
        self.update_db_record(
                    "vpa_analysis",{"_id":self.check_objectid(main_id)}, update_dict)
        insert_text_dict ={
            "name":"name",
            "condition": {"type":"text"}
        }
        update_dict = {
                "text_data": json.dumps(insert_text_dict)}
        self.update_db_record(
                    "vpa_analysis",{"_id":self.check_objectid(main_id)}, update_dict)
        insert_text_list = []
        insert_text = {
            "type" : "text",
            "vpa_id":self.check_objectid(main_id),
            "name" : "Residuals",
            "text" : "Residuals={}".format(self.Residuals),
            "x" : "80%",
            "y" : "95%",
        }
        insert_text_list.append(insert_text)
        self.col_insert_data("vpa_analysis_detail",insert_text_list)

if __name__ == "__main__":
    a = Vpa(None)
    a.add_detail("202104071029108545217538","/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/VPA/env.plot.xls")