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

class PermanovaSingle(ApiBase):
    def __init__(self, bind_object):
        super(PermanovaSingle, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_detail(self, main_id, input_file):
        """
        导入表格数据
        """
        with open(input_file,"r") as infile:
            line = infile.readline()
            fd = line.rstrip().split('\t')
            print(fd)
            title_dict = {
                "Characteristics":"Name",
                "F_Model":"F.Models",
                "p_value":"Pr(>F)",
                "SumsOfSqs":"SumsOfSqs",
                "MeanSqs":"MeanSqs",
                "R2":"R2"
            }
            table_data_insert =[]
            for i in fd:
                if i in title_dict.keys():
                    insert_dict = {
                        "field": i.split('.')[0],
                        "title": title_dict[i],
                        "sort":"false",
                        "filter":"false",
                        "type":"string",
                    }
                    table_data_insert.append(insert_dict)
            permanova_table_dict = {
                "column" : table_data_insert,
                "condition":{'type':"table"}
                }
            update_dict = {
                "table_data": json.dumps(permanova_table_dict)
            }
            self.update_db_record(
                "permanova_single",{"_id":self.check_objectid(main_id)}, update_dict)
            
            while 1:
                insert_table = []
                line = infile.readline()
                if not line:
                    break
                fd1 = line.rstrip().split('\t')
                insert_table_dict = {
                    "type":"table",
                    "permanova_id": self.check_objectid(main_id),
                }
                for i in range(len(fd)):
                    if fd[i] in title_dict.keys():
                        insert_table_dict[fd[i].split('.')[0]] = fd1[i]
                    else:
                        continue
                    print("-----------")
                    print(fd[i])
                    print(fd1[i])
                    print("-"*10)
                insert_table.append(insert_table_dict)
                self.col_insert_data("permanova_single_detail",insert_table)

if __name__ == "__main__":
    a = PermanovaSingle(None)
    a.add_detail("202104071029108545217538","/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/single_permanova/Univariate_PERMANOVA.xls")