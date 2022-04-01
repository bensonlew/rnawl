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


class Bubble(ApiBase):
    def __init__(self, bind_object):
        super(Bubble, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_bubble_detail(self, main_id, target_file):
        with open(target_file, 'r') as r:
            insert_datas = []
            data = r.readlines()
            bubble_data = data[0].strip().split('\t')
            for line in data[1:]:
                line_list = line.strip().split('\t')
                insert_data = {
                    "bubble_id": self.check_objectid(main_id),
                    "type": "bubble",
                }
                for j in range(len(bubble_data)):
                    try:
                        float(line_list[j])
                    except:
                        s_data = line_list[j]
                    else:
                        s_data = float(line_list[j])
                    insert_data[bubble_data[j]] = s_data
                    if 'group' not in bubble_data:
                        insert_data['group'] = ''
                insert_datas.append(insert_data)
        self.col_insert_data("bubble_detail", insert_datas)
        # new_bubble_data = bubble_data
        # if 'name' in bubble_data:
        #     new_bubble_data.remove('name')
        # if 'group' in bubble_data:
        #     new_bubble_data.remove('group')
        bubble_data_insert = {
            "name": "name",
            # "data": new_bubble_data,
            "category": "group",
            "condition": {'type': "bubble"},
        }
        update_dict = {
            "bubble_data": json.dumps(bubble_data_insert)
        }
        self.update_db_record("bubble", {"_id": self.check_objectid(main_id)}, update_dict)


if __name__ == "__main__":
    a = Bubble(None)
    a.add_bubble_detail("5ea7b87917b2bf0f6184ce56", "/mnt/ilustre/users/sanger-dev/workspace/20200507/"
                                                    "Bubble_bubble20200507151717/output/bubble/bubble.xls")
