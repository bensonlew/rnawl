# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
from io import StringIO
import types
import re
import datetime
import os
import json
from api_base import ApiBase

class Yfull(ApiBase):
    def __init__(self, bind_object):
        super(Yfull, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_yfull_detail(self, main_id, text):
        """
        导入结果文件列表
        """
        out_file = text.split('\t')
        insert_file_datas =[]
        # n = 1
        file_info = ""
        for file in out_file:
            if re.search("9000_SNP",file):
                file_info = "样本{}树上已知的位点信息".format(file.split('_')[-3])
            elif re.search("private.xlsx",file):
                file_info = "样本{}私有位点信息".format(file.split('-')[-2])
            elif re.search("private-list",file):
                file_info = "样本{}私有位点信息".format(file.split('-')[-3])
            elif re.search("YFULL_SNP",file):
                file_info = "样本{}}ISOGG收录的位点信息".format(file.split('_')[-3])
            elif re.search("all_sample_result",file):
                file_info = "所有样本分型结果"
            elif re.search("all_sample_dep",file):
                file_info= "所有样本测序深度统计"
            elif re.search("all_sample_MAP_OT",file):
                file_info = "所有样本测序质量统计"
            elif re.search("size.xls",file):
                file_info = "Y染色体bam文件大小统计"
            else:
                file_info = "未知文件"

            insert_data = {
                "type":"table",
                "yfull_id": self.check_objectid(main_id),
                "term": file,
                "result": file_info
            }
            insert_file = []
            insert_file.append(insert_data)
            # n = n + 1
            self.col_insert_data("yfull_detail",insert_file)
        yfull_data_insert = []
        yfull_insert_dict = {
            "field": "term",
            "title": "",
            "sort": "false",
            "filter": "false",
            "type": "string",
        }
        yfull_insert_dict1 = {
            "field": "result",
            "title": "",
            "sort": "false",
            "filter": "false",
            "type": "string",
        }
        yfull_data_insert.append(yfull_insert_dict)
        yfull_data_insert.append(yfull_insert_dict1)
        yfull_table_dict = {
            "column": yfull_data_insert,
            "condition":{'type':'table'}
        }
        update_dict = {
            "yfull_detail_table": json.dumps(yfull_table_dict)
        }
        self.update_db_record(
            "yfull",{"_id": self.check_objectid(main_id)}, update_dict
        )

if __name__ == "__main__":
    a = Yfull(None)
    a.add_yfull_detail("5ea7b87917b2bf0f6184ce88",
                             "all_sample_dep.xls\tall_sample_MAP_OT.xls\tall_sample_result.xlsx\tYG201901399_9000_SNP.xlsx\tYG201901399_YFULL_SNP.xlsx\tYG201901399-private.xlsx\tYG201901399-private-list.xlsx")
