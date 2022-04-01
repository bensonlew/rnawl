# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import json
from api_base import ApiBase

class SingleCorrelation(ApiBase):
    def __init__(self, bind_object):
        super(SingleCorrelation, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_detail(self, main_id, output_file):
        '''
        导入画图表格数据
        '''
        insert_dict = {
            "condition" :{"type":"column"},
            "data":"value",
            "name":"name",
            "category":"category",
        }
        update_dict = {
            "column_data":json.dumps(insert_dict)
        }
        self.update_db_record(
            "single_correlation",{"_id":self.check_objectid(main_id)}, update_dict
        )
        collum_dict = {}
        with open(output_file,"r") as plot:
            while 1:
                line = plot.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                collum_dict[fd[0]] = fd[1]
        sort_value = sorted(zip(collum_dict.itervalues(),collum_dict.iterkeys()),reverse=True)
        sort_key = []
        for i in sort_value:
            sort_key.append(i[1])
        n = 0
        for i in sort_key:
            insert_collumn_data = []
            n += 1
            insert_collumn_dict = {
                "single_co_id":self.check_objectid(main_id),
                "type":"column",
                "name":i,
                "value":float(collum_dict[i]),
                "num": n,
                "category": "positive" if float(collum_dict[i])>0 else "nagetive"
            }
            insert_collumn_data.append(insert_collumn_dict)
            self.col_insert_data("single_correlation_detail",insert_collumn_data)
        
