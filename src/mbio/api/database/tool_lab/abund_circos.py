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

class AbundCircos(ApiBase):
    def __init__(self, bind_object):
        super(AbundCircos, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_chord(self, main_id, taxa_table, percent_table):
        data_list = []
        with open(percent_table,"r") as pt:#获取样本list和orthers
            lines = pt.readlines()
            for line in lines[1:]:
                fd = line.rstrip().split('\t')
                data_list.append(fd[0])
        others = []
        with open(taxa_table,"r") as inf:
            line = inf.readline()
            feators1 = []
            feator_data1 = {}
            feators = []
            feator_data = {}
            fd = line.rstrip().split('\t')
            for i in fd[1:]:
                feators1.append(i)
                feator_data1[i] = [0]*len(fd[1:])
                others = [0]*len(fd[1:])
                feator_data[i] = []
            n = 1
            while 1:
                line1 = inf.readline()
                if not line1:
                    break
                fd1 = line1.rstrip().split('\t')
                if fd1[0] in data_list:
                    n+=1
                    feators.append(fd1[0])
                    feator_data[fd1[0]] = []
                    feator_data1[fd1[0]] = []
                    data = fd1[1:]
                    for i in range(len(data)):
                        feator_data[feators1[i]].append(float(data[i]))
                        feator_data1[fd1[0]].append(float(data[i]))
                else:
                    data = fd1[1:]                    
                    for i in range(len(data)):
                        others[i] += float(data[i])  
            feators.append("others")                
            feators.extend(feators1)
            temp_others = [0]*n
            temp_others.extend(others)
            feator_data["others"] = temp_others
            for yf in range(len(others)):
                feator_data[fd[yf+1]].append(others[yf])
            for j in feators:
                if j in fd[1:]:
                    feator_data[j].extend(feator_data1[j])
                    continue
                if j == "others":
                    continue
                feator_data[j].extend([0]*n)
                feator_data[j].extend(feator_data1[j])
            insert_head_dict = {
                "circos_id":self.check_objectid(main_id),
                "type":"chord",
                "index":1,
                "data":feators
            }
            insert_head_list = [insert_head_dict]
            self.col_insert_data("abund_circos_chord",insert_head_list)
            num = 2
            for n in feators:
                insert_dict = {
                    "circos_id":self.check_objectid(main_id),
                    "type":"chord",
                    "index":num,
                    "data":feator_data[n]
                }
                insert_chord_list  = [insert_dict]
                self.col_insert_data("abund_circos_chord",insert_chord_list)
                num+=1
        insert_dict = {
            "data":"data",
            "condition":{"type":"chord"}
        }
        update_dict = {
            "chord_data":json.dumps(insert_dict)
        }
        self.update_db_record(
            "abund_circos",{"_id":self.check_objectid(main_id)}, update_dict
        )
    
    def add_arc(self,main_id,percent_table,group_table):
        # data2 = []
        with open(percent_table,"r") as tpt:
            line = tpt.readline()
            fd = line.rstrip().split('\t')
            group = fd[1:]
            group_index = {}
            n = 0
            for i in fd[1:]:
        
                group_index[n] = i
                n += 1
            while 1:
                line = tpt.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                data = fd[1:]
                for i in range(len(data)):
                    dict_temp = {
                        "circos_id":self.check_objectid(main_id),
                        "type": "chord_arc",
                        "name" : fd[0],
                        "group" : group_index[i],
                        "value" : float(data[i]),
                        "index" : 1
                    }
                    insert_list = [dict_temp]
                    self.col_insert_data("abund_circos_arc",insert_list)
                    # data2.append(dict_temp)

        with open(group_table,"r") as tpt:
            line = tpt.readline()
            fd = line.rstrip().split('\t')
            group = fd[1:]
            group_index = {}
            n = 0
            for i in fd[1:]:
        
                group_index[n] = i
                n += 1
            while 1:
                line = tpt.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                data = fd[1:]
                for i in range(len(data)):
                    dict_temp = {
                        "circos_id":self.check_objectid(main_id),
                        "type": "chord_arc",
                        "name" : group_index[i],
                        "group" : fd[0],
                        "value" : float(data[i]),
                        "index" : 2
                    }
                    insert_list = [dict_temp]
                    self.col_insert_data("abund_circos_arc",insert_list)
                    # data2.append(dict_temp)
        insert_dict = {
            "name":"name",
            "index":"index",
            "value":"value",
            "group":"group",
            "condition":{"type":"chord_arc"}
        }
        update_dict = {
            "chord_arc_data":json.dumps(insert_dict)
        }
        self.update_db_record(
            "abund_circos",{"_id":self.check_objectid(main_id)}, update_dict
        )
            

            

            