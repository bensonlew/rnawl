# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'wuqin'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import types
import re
import datetime
import os
import json
from api_base import ApiBase


class HistogramApi(ApiBase):
    """频率直方图导表"""
    def __init__(self, bind_object):
        super(HistogramApi, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_histogram_detail(self, main_id, histogram_table):
        """
        频率直方图详情表导入
        """
        # histogram_table = os.path.join(output_dir + 'histogram/histogram.csv')

        with open(histogram_table, 'r') as r:
            insert_datas = []
            data = r.readlines()
            histogram_data = ['name', 'value']
            for line in data:
                line_list = line.strip().split('\t')
                insert_data = {
                    "histogram_id": self.check_objectid(main_id),
                    "type": "column",
                }
                for j in range(len(histogram_data)):
                    try:
                        float(line_list[j])
                    except:
                        h_data = line_list[j]
                    else:
                        h_data = float(line_list[j])
                    insert_data[histogram_data[j]] = h_data
                    insert_data['category'] = ''
                insert_datas.append(insert_data)
        self.col_insert_data("histogram_detail", insert_datas)
        histogram_data_insert = {
            "name": "name",
            "data": 'value',
            "category": "category",
            "condition": {'type': "column"}
        }
        update_dict = {
            "column_data": json.dumps(histogram_data_insert)
        }
        self.update_db_record("histogram", {"_id": self.check_objectid(main_id)}, update_dict)


if __name__ == "__main__":
    a = HistogramApi(None)
    a.add_histogram_detail("5d6c843f17b2bf4bad32f53b", "/mnt/ilustre/users/sanger-dev/workspace/"
                                                       "20200420/Single_histogram20200420123953/Histogram/output")
