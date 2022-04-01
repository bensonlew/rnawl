# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import gzip
import re
import datetime
import os
import json
from api_base import ApiBase

class VcfLddecay(ApiBase):
    def __init__(self, bind_object):
        super(VcfLddecay, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_detail(self, main_id, picture_path):
        """
        基于ld-decay分析
        """
        insert_dict = {
            "picture":picture_path
         }
        self.update_db_record(
            "vcf_lddecay",{"_id":self.check_objectid(main_id)},insert_dict
        )
    
    def vcf_ld_detail(self, main_id, file_path):
        ld_id = self.check_objectid(main_id)
        with open(file_path)as r:
            header = []
            dict1 = {}
            data_list1 = []
            data_list2 = []
            data_list3 = []
            lines = r.readlines()
            for i in lines[1:]:
                list1 = []
                group_name = i.strip().split("\t")[0]
                header.append(group_name)
                r2_8_index = 0  # 第一列的数字刚开始赋值为0。
                r2_8_data = 3  # 这三个字都赋值为-3，保证肯定不会取到这三个值。
                r2_1_index = 0
                r2_1_data = 3
                r2_5_index = 0
                r2_5_data = 3
                max1 = 0
                half_life = 0
                with gzip.open(i.strip().split("\t")[1], "r") as f:
                    ld_lines = f.readlines()
                    for k in ld_lines[1:]:
                        if float(k.split("\t")[1]) > max1:
                            max1 = float(k.split("\t")[1])
                            half_life = round(max1 / float(2), 4)
                    for j in ld_lines[1:]:
                        if float(j.split("\t")[1])-float(0.8) >= 0:
                            if float(r2_8_data)-float(0.8) > float(j.split("\t")[1]) - float(0.8):
                                r2_8_data = j.split("\t")[1]
                                r2_8_index = j.split("\t")[0]
                        if float(j.split("\t")[1]) - float(0.1) >= 0:
                            if float(r2_1_data) - float(0.1) > float(j.split("\t")[1]) - float(0.1):
                                r2_1_data = j.split("\t")[1]
                                r2_1_index = j.split("\t")[0]
                        if float(j.split("\t")[1]) - half_life >= 0:
                            if float(r2_5_data) - half_life > float(j.split("\t")[1]) - half_life:
                                r2_5_data = j.split("\t")[1]
                                r2_5_index = j.split("\t")[0]
                    if r2_8_index == 0:
                        r2_8_index = "/ "
                    list1.append(r2_8_index)
                    if r2_1_index == 0:
                        r2_1_index = "/ "
                    list1.append(r2_1_index)
                    if r2_5_index == 0:
                        r2_5_index = "/ "
                    list1.append(r2_5_index)
                dict1[i.strip().split("\t")[0]] = list1
            table_r28 = {}
            table_r21 = {}
            table_r25 = {}
            for x in header:
                table_r28[x] = dict1[x][0]
                table_r21[x] = dict1[x][1]
                table_r25[x] = dict1[x][2]
            # for n in ["r2_8", "r2_1", "r2_5"]:
            for n in ["R2 > 0.8", "R2 > 0.1", "R2 half-life"]:
                insert_data = {"pops": n, "ld_id": self.check_objectid(main_id),"type":"table"}
                if n == "R2 > 0.8":
                    for key1 in table_r28.keys():
                        insert_data[key1] = table_r28[key1]
                    data_list1.append(insert_data)
                    if len(data_list1) == 0:
                        self.bind_object.logger.info("r2_8{}文件为空！".format(file_path))
                    else:
                        self.col_insert_data("vcf_lddecay_detail", data_list1)
                elif n == "R2 > 0.1":
                    for key2 in table_r21.keys():
                        insert_data[key2] = table_r21[key2]
                    data_list2.append(insert_data)
                    if len(data_list2) == 0:
                        self.bind_object.logger.info("r2_1{}文件为空！".format(file_path))
                    else:
                        self.col_insert_data("vcf_lddecay_detail", data_list2)
                elif n == "R2 half-life":
                    for key3 in table_r25.keys():
                        insert_data[key3] = table_r25[key3]
                    data_list3.append(insert_data)
                    if len(data_list3) == 0:
                        self.bind_object.logger.info("r2_1{}文件为空！".format(file_path))
                    else:
                        self.col_insert_data("vcf_lddecay_detail", data_list3)
            table_data_insert = [{
                "field":"pops",
                "title":"POPs",
                "sort":"false",
                "filter":"false",
                "type":"string"
            }]
            for i in header:
                insert_dict = {
                    "field":i,
                    "title":i,
                    "sort":"false",
                    "filter":"false",
                    "type":"string"
                }
                table_data_insert.append(insert_dict)
            ld_table_dict = {
                "column" : table_data_insert,
                "condition": {"type":"table"}
            }
            update_dict = {
                "table_data":json.dumps(ld_table_dict)
            }
            self.update_db_record("vcf_lddecay", {"_id": ld_id},update_dict )