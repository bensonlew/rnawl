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
import pandas as pd
import json
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from bson.objectid import ObjectId


class AbundanceBubble(ApiBase):
    def __init__(self, bind_object):
        super(AbundanceBubble, self).__init__(bind_object)
        self._project_type = "tool_lab"

    def add_abundance_bubble_detail(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None, abundance_table=None):

        if main_id is None:
            name = "abundance_bubble" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='abundance_bubble',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('abundance_bubble', [main_info])
        else:
            main_id = ObjectId(main_id)

        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        f = pd.read_table(abundance_table, sep="\t")
        max_num = f.ix[1:, 1:].max(axis=1).max()
        min_num = f.ix[1:, 1:].min(axis=1).min()
        median_num = f.ix[1:, 1:].median(axis=0).median()
        customSize = str([min_num, median_num, max_num])

        data_list_abundance = []
        with open(abundance_table) as f1:
            lines = f1.readlines()
            group_list = lines[0].strip().split('\t')[1:]
            for line in lines[1:]:
                item = line.strip().split("\t")
                for i in range(len(group_list)):
                    insert_data = {
                        "main_id": main_id,
                        "name": item[0],
                        "x": group_list[i],
                        "y": item[0],
                        "size": round(float(item[i+1]), 2),
                        "color": item[0],
                        "group": "",
                        "type": "bubble",
                    }
                    data_list_abundance.append(insert_data)
        self.create_db_table("abundance_bubble_detail", data_list_abundance)
        update_dict = json.dumps({
            "name": "name",
            "fdr": "color",
            "category": "group",
            "condition": {'type': "bubble"}
        })

        self.update_db_record("abundance_bubble", main_id, main_id=main_id, bubble_data=update_dict, customSize=customSize, status="end")


    def add_abundance_bubble_detail2(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None, abundance_table=None, taxonomy_table=None):

        if main_id is None:
            name = "abundance_bubble" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='abundance_bubble',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('abundance_bubble', [main_info])
        else:
            main_id = ObjectId(main_id)

        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        f = pd.read_table(abundance_table, sep="\t")
        max_num = f.ix[1:, 1:].max(axis=1).max()
        min_num = f.ix[1:, 1:].min(axis=1).min()
        median_num = f.ix[1:, 1:].median(axis=0).median()
        customSize = str([min_num, median_num, max_num])

        taxonomy_dict = {}
        data_list_abundance = []
        with open(abundance_table) as f1, open(taxonomy_table) as f2:

            lines = f1.readlines()
            group_list = lines[0].strip().split('\t')[1:]
            rows = f2.readlines()
            for row in rows[1:]:
                item = row.strip().split('\t')
                if item[0] not in taxonomy_dict.keys():
                    taxonomy_dict[item[0]] = item[1]
            for line in lines[1:]:
                item = line.strip().split("\t")
                if item[0] in taxonomy_dict.keys():
                    taxonomy = taxonomy_dict[item[0]]

                for i in range(len(group_list)):
                    insert_data = {
                        "main_id": main_id,
                        "name": item[0],
                        "x": group_list[i],
                        "y": item[0],
                        "size": round(float(item[i+1]), 2),
                        "color": taxonomy,
                        "group": "",  # 20210916 group存空
                        "type": "bubble",
                    }
                    data_list_abundance.append(insert_data)
        self.create_db_table("abundance_bubble_detail", data_list_abundance)

        update_dict = json.dumps({
                "name": "name",
                "fdr": "color",
                "category": "group",
                "condition": {'type': "bubble"}
            })

        self.update_db_record("abundance_bubble", main_id, main_id=main_id, bubble_data=update_dict, customSize=customSize, status="end")

