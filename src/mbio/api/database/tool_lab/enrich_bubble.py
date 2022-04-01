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

class EnrichBubble(ApiBase):
    def __init__(self, bind_object):
        super(EnrichBubble, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_detail(self, main_id, target_dir, major):
        groups = []
        for root, dir_names, files in os.walk(target_dir):
            x = 0
            scale = {}
            for file in files:    
                with open(os.path.join(root,file),"r") as tb:
                    group_name = file[:-11]
                    insert_datas = []
                    data = tb.readlines()
                    bubble_data = data[0].strip().split('\t')
                    for line in data[1:]:
                        x = x +1
                        line_list = line.strip().split('\t')
                        insert_data = {
                            "name" : "x{}".format(x),
                            "bubble_id": self.check_objectid(main_id),
                            "type": "bubble",
                        }
                        if scale.has_key("min"):
                            if scale["min"]>int(line_list[2]):
                                scale["min"] = int(line_list[2])
                        else:
                            scale["min"] = int(line_list[2])
                        if scale.has_key("max"):
                            if scale["max"]<int(line_list[2]):
                                scale["max"] = int(line_list[2])
                        else:
                            scale["max"] = int(line_list[2])
                        for j in range(len(bubble_data)):
                            try:
                                float(line_list[j])
                            except:
                                s_data = line_list[j]
                            else:
                                s_data = float(line_list[j])
                            insert_data[bubble_data[j].lower()] = s_data
                            insert_data["group"] = group_name
                            if group_name != major:
                                insert_data['num'] = 0
                        insert_datas.append(insert_data)
                    if len(insert_datas) > 0:
                        groups.append(group_name)
                        self.col_insert_data("enrich_bubble_detail", insert_datas)
        bubble_data_insert = {
            "name": "name",
                # "data": new_bubble_data,
            "category": "group" ,
            "fdr":"color",
            "condition": {'type': "bubble"},
        }
        update_dict = {
            "bubble_data": json.dumps(bubble_data_insert)
            }
        self.update_db_record("enrich_bubble", {"_id": self.check_objectid(main_id)}, update_dict)
        group_dict = {
            "group": groups
        }
        self.update_db_record("enrich_bubble", {"_id": self.check_objectid(main_id)}, group_dict)
        if len(scale.keys()) == 0:
            scale['min'] = 0
            scale['max'] = 0
        tem_scale = (scale["max"] - scale["min"])/3
        scale1 = scale["max"] - tem_scale
        scale2 = scale['min'] + tem_scale
        scale_data_insert = [scale['min'],scale2,scale1,scale["max"]]
        update_dict = {
            "customizeSize": json.dumps(scale_data_insert)
            }
        self.update_db_record("enrich_bubble", {"_id": self.check_objectid(main_id)}, update_dict)

    def add_tooltip(self,main_id,file, method):
        p_name = ""
        term_name = ""
        with open(file,"r") as finfo:
            line = finfo.readline()
            fd =line.rstrip().split('\t')
            p_name = fd[0].capitalize()
            if method == "go":
                term_name = "GO Term"
            elif method == "kegg":
                term_name = "Pathway"
        
        group_data_insert = {
            "x": "Rich Factor",
            "y": term_name ,
            "size": "Num" ,
            "fdr": p_name,
        }
        update_dict = {
            "tooltip_name": json.dumps(group_data_insert)
            }
        self.update_db_record("enrich_bubble", {"_id": self.check_objectid(main_id)}, update_dict)