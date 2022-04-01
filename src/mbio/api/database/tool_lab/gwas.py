# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

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
from api_base import ApiBase


class Gwas(ApiBase):
    def __init__(self, bind_object):
        super(Gwas, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_gwas_main(self, result, pic_data, trait, main_id=None, params=None, project_sn='gwas', task_id='gwas'):
        # add main table info
        if main_id is None:
            name = "GWAS" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                desc='GWAS',
                params=params,
                version="v1",
                status="start",
            )
            main_id = self.create_db_table('sg_gwas', [main_info])
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)

        # detail_data required to display table
        detail_data = {
            'column': [],
            'condition': {},
        }
        page_col = ['SNP', 'Chromosome', 'Base Pair', 'Effect', 'P Value']
        data_col = ['snp', 'chrom', 'bp', 'effect', 'pvalue']
        length = len(data_col)
        for i in range(length):
            if i in [0, 1]:
                d_type = 'string'
            elif i in [2]:
                d_type = 'int'
            else:
                d_type = 'float'
            col_detail = {
                "filter": "false",
                "field": data_col[i],
                "title": page_col[i],
                "type": d_type,
                "sort": "false",
            }
            detail_data['column'].append(col_detail)

        trait_d = {'data': trait}
        query_dict = {
            '_id': main_id
        }
        update_dict = {
            'status': 'end',
            'main_id': main_id,
            'detail_data': json.dumps(detail_data, sort_keys=True, separators=(',', ':')),
            'pic_data': pic_data,
            'traits': json.dumps(trait_d, sort_keys=True, separators=(',', ':')),
        }
        for each in result:
            if os.path.exists(each):
                self.add_gwas_detail(main_id=main_id, data=each)
        self.update_db_record('sg_gwas', query_dict=query_dict, update_dict=update_dict)
        return main_id

    def add_gwas_detail(self, main_id, data):
        trait = os.path.basename(data).split('_gwas_loci')[0].split('.')
        if len(trait) > 1:
            trait = '_'.join(trait)
        else:
            trait = trait[0]
        data_df = pd.read_table(data, header=0, dtype={'snp': 'str', 'chrom': 'str'})
        data_df["gwas_id"] = main_id
        data_df['trait'] = trait
        gwas_list = data_df.to_dict('records')
        try:
            self.col_insert_data('sg_gwas_detail', gwas_list)
        except Exception as e:
            self.bind_object.logger.info("导入sg_gwas_detail:%s信息出错:%s" % (main_id, e,))
        else:
            self.bind_object.logger.info("导入sg_gwas_detail:%s成功" % (main_id,))