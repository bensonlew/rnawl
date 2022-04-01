# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20200812
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

class Ssgsea(ApiBase):
    def __init__(self, bind_object):
        super(Ssgsea, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_ssgsea(self, main_id, es_matrix, gmx, params=None, project_sn='medical_transcriptome', task_id='medical_transcriptome'):
        """
        dump exp data into database
        :param exp_matrix: expression matrix path or an express matrix in pandas DataFrame format.
        :param exp_level: str, transcript or gene
        :param exp_type: str, usually is tpm or fpkm or count. default tpm
        :param group_dict: ordered dict of group info
        :param quant_method: the method to be used to quant expression.
        :param project_sn: project id
        :param task_id: task id
        :param main_id: 主表id，如果提供，则本函数不创建主表
        :param group_id: 包含分组方案信息的主表id
        :param add_distribution: 是否添加表达分布信息到数据库
        :param params: parameters dict for expression quant.
        :return: main table id
        :version:version info
        """
        # add main table info
        if main_id is None:
            name = "SSgsea" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='SSgsea',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_geneset_ssgsea', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        df = pd.read_table(es_matrix, header=0, sep='\t')
        gsva_columns = df.columns
        columns_list = list()
        columns_list.append({'field': 'geneset', 'filter': False, 'sort':False, 'title': 'geneset', 'type': 'string'})
        for i in gsva_columns[1:]:
            columns_list.append({'field': i, 'filter': False, 'sort': False, 'title': i, 'type': 'float'})
        data_columns = {'column': columns_list, 'condition': {}}
        column_data = json.dumps(data_columns)
        df['ssgsea_id'] =main_id
        detail = df.to_dict('r')
        self.create_db_table('sg_geneset_ssgsea_detail', detail)

        main_collection = self.db['sg_geneset_ssgsea']
        main_collection.update({"_id": main_id},
                               {"$set": {"status": "end", 'column_data': column_data}})
        return main_id
