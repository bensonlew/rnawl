# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20210303
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class RrnaStat(ApiBase):
    def __init__(self, bind_object):
        super(RrnaStat, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_rrna_stat(self, main_id, rrna_stat_table):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(rrna_stat_table, header=0, sep='\t')
        df.columns=["sample","ratio"]
        df["rrna_stat_id"] = main_id
        columns_list = list()
        columns_list.append({'field': 'sample', 'filter': False, 'sort': False, 'title': 'Sample', 'type': 'str'})
        columns_list.append({'field': 'ratio', 'filter': False, 'sort': False, 'title': 'rRNA(%)', 'type': 'str'})
        data_columns = {'column': columns_list, 'condition': {}}
        columns_data = json.dumps(data_columns)
        self.create_db_table('sg_rrna_stat_detail', df.to_dict('r'))
        self.update_db_record('sg_rrna_stat', main_id, main_id=main_id, column_data_detail=columns_data)

