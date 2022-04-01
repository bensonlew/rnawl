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

class Sankey(ApiBase):
    def __init__(self, bind_object):
        super(Sankey, self).__init__(bind_object)
        self._db_name = "sanger_tool_lab"

    def add_detail(self, main_id, input_file):
        """
        导入桑基图数据
        """
        with open(input_file,"r") as infile:
            insert_dict = {
                "condition":{"type":"sankey_new"},
                "source":"source",
                "target":"target",
                "value":"value"
            }
            update_dict = {
                "sankey_new_data": json.dumps(insert_dict)
            }
            self.update_db_record(
                "sankey",{"_id":self.check_objectid(main_id)}, update_dict)
            while 1:
                insert_sankey = []
                line = infile.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                source = fd[0]
                target = fd[1]
                value = float(fd[2])
                insert_table_dict ={
                    "type":"sankey_new",
                    "sankey_id":self.check_objectid(main_id),
                    "source":source,
                    "target":target,
                    "value":value
                }
                insert_sankey.append(insert_table_dict)
                self.col_insert_data("sankey_detail",insert_sankey)
        


                
