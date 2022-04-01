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

class VcfTree(ApiBase):
    def __init__(self, bind_object):
        super(VcfTree, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_tree_detail(self, main_id, file_path, group_file):
        with open(file_path, "r") as r:
            data = r.readlines()
            tree_str = data[0].strip()
        group = {}
        with open(group_file, "r") as gf:
            while 1:
                line = gf.readline()
                if not line:
                    break
                fd = line.rstrip().split("\t")
                if not group.has_key(fd[1]):
                    group[fd[1]] = [fd[0]]
                else:
                    group[fd[1]].append(fd[0])
        insert_group_list = []
        for i in sorted(group.keys()):
            insert_group_dict = {
                "groupname":i,
                "value": group[i]
            }
            insert_group_list.append(insert_group_dict)
        insert_tree_dict = {
            "tree_id":self.check_objectid(main_id),
            "name":"h_tree",
            "data": tree_str,
            "direction":"h",
            "type":"tree",
            "group":insert_group_list
        }
        insert_tree_dict1 = {
            "tree_id":self.check_objectid(main_id),
            "name":"c_tree",
            "data": tree_str,
            "direction":"circle",
            "type":"tree",
            "group":insert_group_list
        }
        insert_tree_list=[insert_tree_dict,insert_tree_dict1]
        self.col_insert_data("vcf_tree_detail",insert_tree_list)
        insert_dict ={
            "name":"name",
            "condition":{"type":"tree"}
        }
        update_dict = {
            "tree_data":json.dumps(insert_dict)
        }
        self.update_db_record(
            "vcf_tree",{"_id":self.check_objectid(main_id)}, update_dict
        )
        insert_dict1= {
            "data":["circle","h"],
        }
        update_dict1 = {
            "direction_data" : json.dumps(insert_dict1)
        }
        self.update_db_record(
            "vcf_tree",{"_id":self.check_objectid(main_id)}, update_dict1
        )
            
                