# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'fangyifei'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
import types
import re
import datetime
import os
from bson.objectid import ObjectId
import json
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class MultigroupBar(ApiBase):
    def __init__(self, bind_object):
        super(MultigroupBar, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_multigroup_bar_detail(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None, result_file=None, combination=None):
        if main_id is None:
            name = "multigroup_bar" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='multigroup_bar',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('multigroup_bar', [main_info])
        else:
            main_id = ObjectId(main_id)
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        data_list_bar = []
        data_list_ishape = []
        with open(result_file) as f:
            lines = f.readlines()
            lines.reverse()  # 20210918 修改为倒序
            for line in lines[:-1]:
                item = line.strip().split("\t")
                if combination != "none":
                    insert_data = {
                        "main_id": main_id,
                        "name": item[0].replace('"', ''),
                        "category": item[3],
                        "value": float(item[1]),
                        "type": "column",
                    }
                    data_list_bar.append(insert_data)

                    insert_data2 = {
                        "main_id": main_id,
                        "name": item[0].replace('"', ''),
                        "mean": float(item[1]),
                        "sd_high": float(item[2]),
                        "sd_low": float(item[2]),
                        "category": item[3],
                        "type": "ishape",
                    }
                    data_list_ishape.append(insert_data2)
                else:
                    insert_data = {
                        "main_id": main_id,
                        "name": item[0].replace('"', ''),
                        "category": item[2],
                        "value": float(item[1]),
                        "type": "column",
                    }
                    data_list_bar.append(insert_data)

        self.create_db_table("multigroup_bar_detail", data_list_bar)
        if combination != "none":
            self.create_db_table("multigroup_bar_ishape", data_list_ishape)

        update_dict1 = json.dumps({
                "name": "name",
                "data": "value",
                "category": "category",
                "condition": {'type': "column"}
            })

        update_dict2 = json.dumps({
                "name": "name",
                "data": ["mean", "sd_high", "sd_low"],
                "group": "category",
                "condition": {'type': "ishape"}
            })

        self.update_db_record("multigroup_bar", main_id, main_id=main_id, column_data=update_dict1, ishape_data=update_dict2, status="end")
        # if combination != "none":  20210918 详情表结构必须保持固定格式
            # self.update_db_record("multigroup_bar", main_id, main_id=main_id, ishape_data=update_dict2, status="end")






