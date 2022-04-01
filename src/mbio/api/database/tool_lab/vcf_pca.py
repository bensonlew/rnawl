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

class VcfPca(ApiBase):
    def __init__(self,bind_object):
        super(VcfPca, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"
    
    def add_detail(self, main_id, val_file, vec_file, group_file):
        """
        导入PCA图数据
        """
        sample_group = {}
        axis = []
        with open(group_file,"r") as gp:
            while 1:
                line = gp.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                sample_group[fd[0]]=fd[1]
        with open(vec_file, "r") as vec:
            while 1:
                insert_scatter = []
                insert_table = []
                line = vec.readline()
                if not line:
                    break
                fd = line.rstrip().split(' ')
                data = fd[2:]
                insert_scatter_dict = {
                    "type" : "scatter",
                    "pca_id": self.check_objectid(main_id),
                    "group": sample_group[fd[0]],
                    "name": fd[0],
                }
                insert_table_dict = {
                    "type":"table",
                    "pca_id": self.check_objectid(main_id),
                    "sample_id":fd[0]
                }
                for i in range(len(data)):
                    insert_scatter_dict["PC{}".format(i+1)] = float(fd[i+2])
                    insert_table_dict["PC{}".format(i+1)] = fd[i+2]
                insert_scatter.append(insert_scatter_dict)
                insert_table.append(insert_table_dict)
                self.col_insert_data("vcf_pca_scatter",insert_scatter)
                self.col_insert_data("vcf_pca_table",insert_table)
        for i in range(1,21):
            axis.append("PC{}".format(i))
        insert_dict = {
            "category":"group",
            "data":axis,
            "name":"name",
            "condition":{"type":"scatter"},
        }
        update_dict = {
            "scatter_data": json.dumps(insert_dict)
        }
        self.update_db_record(
            "vcf_pca",{"_id":self.check_objectid(main_id)},update_dict
        )
        table_data_insert= [{
                "field": "sample_id",
                "title": "sample_ID",
                "sort":"false",
                "filter":"false",
                "type":"string",
            }]
        for i in axis:
            tem_insert_dict = {
                "field": i,
                "title": i,
                "sort":"false",
                "filter":"false",
                "type":"string",
            }
            table_data_insert.append(tem_insert_dict)
        insert_dict1 ={
            "column":table_data_insert,
            "condition":{"type":"table"}
        }
        update_dict1 = {
            "table_data": json.dumps(insert_dict1)
        }
        self.update_db_record(
            "vcf_pca",{"_id":self.check_objectid(main_id)}, update_dict1
        )
        with open(val_file,"r") as vl:
            n = 1
            while 1:
                insert_table = []
                line = vl.readline()
                if not line:
                    break
                if n > 20:
                    break
                value = line.rstrip()
                insert_table_dict ={
                    "type":"table",
                    "pca_id": self.check_objectid(main_id),
                    "name": "PC{}".format(n),
                    "eigenval":value
                }
                n+=1
                insert_table.append(insert_table_dict)
                self.col_insert_data("vcf_pca_eigenval",insert_table)
        insert_dict2 = {
            "column":[{
                "field": "name",
                "title": "name",
                "sort":"false",
                "filter":"false",
                "type":"string",
            },{
                "field": "eigenval",
                "title": "eigenval",
                "sort":"false",
                "filter":"false",
                "type":"string",
            }],
            "condition":{"type":"table"}
        }
        update_dict2 = {
            "eigenval_data": json.dumps(insert_dict2)
        }
        self.update_db_record(
            "vcf_pca",{"_id":self.check_objectid(main_id)}, update_dict2
        )

        



                
                




