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


class Scatter(ApiBase):
    def __init__(self, bind_object):
        super(Scatter, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def add_scatter_detail(self, main_id, target_file):
        # target_file = os.path.join(output_dir, 'scatter.xls')
        with open(target_file, 'r') as r:
            insert_datas = []
            data = r.readlines()
            scatter_data = data[0].strip().split('\t')
            for line in data[1:]:
                line_list = line.strip().split('\t')
                insert_data = {
                    "scatter_id": self.check_objectid(main_id),
                    "type": "scatter",
                }
                for j in range(len(scatter_data)):
                    try:
                        float(line_list[j])
                    except:
                        s_data = line_list[j]
                    else:
                        s_data = float(line_list[j])
                    insert_data[scatter_data[j]] = s_data
                    if 'group' not in scatter_data:
                        insert_data['group'] = ''
                insert_datas.append(insert_data)
        self.col_insert_data("scatter_detail", insert_datas)
        new_scatter_data = scatter_data
        if 'name' in scatter_data:
            new_scatter_data.remove('name')
        if 'group' in scatter_data:
            new_scatter_data.remove('group')
        scatter_data_insert = {
            "name": "name",
            "data": new_scatter_data,
            "category": "group",
            "condition": {'type': "scatter"},
        }
        update_dict = {
            "scatter_data": json.dumps(scatter_data_insert)
        }
        self.update_db_record("scatter", {"_id": self.check_objectid(main_id)}, update_dict)


if __name__ == "__main__":
    a = Scatter(None)
    a.add_scatter_detail("5ea7b87917b2bf0f6184ce56", "/mnt/ilustre/users/sanger-dev/workspace/20200428/Scatter_tsg_"
                                                     "3421_0428130041915383_3852/output/scatter/scatter.xls")
