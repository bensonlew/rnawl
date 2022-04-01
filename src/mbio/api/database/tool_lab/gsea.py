# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:20190604
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class Gsea(ApiBase):
    def __init__(self, bind_object):
        super(Gsea, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def run_webroot(self, main_id, gsea_dir, result_path, term=None):
        task_id = self.bind_object.sheet.id
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或者其对应的字符串！", code="53703606")

        self.bind_object.logger.info("开始导入gsea结果")
        for path in ['gsea_report.xls', 'all_exp.detail', 'all_sets.detail']:
            if os.path.exists(os.path.join(gsea_dir, path)):
                pass
            else:
                self.bind_object.set_error("%s 不存在", variables=(os.path.join(gsea_dir, path)), code="53703607")
        gsea_id = ObjectId(main_id)
        gsea_columns = {'column': [{'field': 'gene_set_name', 'filter': False, 'sort': False, 'title': 'gene_set_name',
                                    'type': 'string'},
                                   {'field': 'description', 'filter': False, 'sort': False,
                                    'title': 'description',
                                    'type': 'string'},
                                   {'field': 'group', 'filter': False, 'sort': False, 'title': 'group',
                                    'type': 'string'},
                                   {'field': 'size', 'filter': False, 'sort': False, 'title': 'size',
                                    'type': 'int'},
                                   {'field': 'es', 'filter': False, 'sort': False, 'title': 'es',
                                    'type': 'float'},
                                   {'field': 'nes', 'filter': False, 'sort': False, 'title': 'nes',
                                    'type': 'float'},
                                   {'field': 'p_value', 'filter': False, 'sort': False, 'title': 'p_value',
                                    'type': 'float'},
                                   {'field': 'fdr_q_value', 'filter': False, 'sort': False, 'title': 'fdr_q_value',
                                    'type': 'float'},
                                   {'field': 'fwer_p_value', 'filter': False, 'sort': False, 'title': 'fwer_p_value',
                                    'type': 'float'},
                                   {'field': 'rank_at_max', 'filter': False, 'sort': False, 'title': 'rank_at_max',
                                    'type': 'int'},
                                   {'field': 'leading_edge_num', 'filter': False, 'sort': False, 'title': 'leading_edge_num',
                                    'type': 'int'},
                                   {'field': 'leading_edge_genes', 'filter': False, 'sort': False, 'title': 'leading_edge_genes',
                                    'type': 'string'},
                                   {'field': 'path', 'filter': False, 'sort': False,
                                    'title': 'path',
                                    'type': 'string'}], 'condition': {}}
        column_data = json.dumps(gsea_columns)
        self.add_gsea_stat(gsea_id, os.path.join(gsea_dir, 'gsea_report.xls'), result_path)
        select_field = [{'field': 'fdr_q_value', 'title': 'Padjust'}, {'field': 'p_value', 'title': 'P value'}]
        select_field_data = json.dumps(select_field)
        select_operate = [{"field": ">", "title": "大于"}, {"field": "<", "title": "小于"},{"field": "=", "title": "等于"},
                          {"field": "<=", "title": "小于等于"},{"field": ">=", "title": "大于等于"}]
        select_operate_data = json.dumps(select_operate)
        # self.add_gsea_geneset(gsea_id, os.path.join(gsea_dir, 'all_sets.detail'))
        # self.add_gsea_exp(gsea_id, os.path.join(gsea_dir, 'all_exp.detail'))
        main_collection = self.db['sg_geneset_gsea']
        main_collection.update({"_id": main_id},
                               {"$set": {"gsea_sets": term , "status": "end", 'column_data': column_data,
                                         'select_field': select_field_data,
                                         'select_operate': select_operate_data}})


    @report_check
    def add_gsea_stat(self, gsea_id, gsea_stat, result_path):
        gsea_stat_df = pd.read_table(gsea_stat, header=0)
        gsea_stat_df = gsea_stat_df.fillna("")
        gsea_stat_df["gsea_id"] = gsea_id
        gsea_stat_df['path'] = result_path
        row_dict_list = gsea_stat_df.to_dict('records')
        print row_dict_list
        main_collection = self.db['sg_geneset_gsea_stat']
        self.create_db_table('sg_geneset_gsea_stat', row_dict_list)
        # try:
        #     self.create_db_table('sg_geneset_gsea_stat', row_dict_list)
        # except Exception as e:
        #     self.bind_object.set_error("导入main: %s出错!" , variables=(gsea_stat), code="53703608")
        # else:
        #     self.bind_object.logger.info("导入gsea_stat：%s出错!" % (gsea_stat))

    @report_check
    def add_gsea_geneset(self, gsea_id, gsea_geneset):
        gsea_geneset_df = pd.read_table(gsea_geneset, header=0)
        gsea_geneset_df['gsea_id'] = gsea_id
        row_dict_list = gsea_geneset_df.to_dict('records')
        main_collection = self.db['sg_geneset_gsea_geneset']

        try:
            self.create_db_table('sg_geneset_gsea_geneset', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=(gsea_geneset), code="53703609")
        else:
            self.bind_object.logger.info("导入gsea_geneset：%s出错!" % (gsea_geneset))

    @report_check
    def add_gsea_exp(self, gsea_id, gsea_exp):
        gsea_exp_df = pd.read_table(gsea_exp, header=0)
        gsea_exp_df['gsea_id'] = gsea_id
        row_dict_list = gsea_exp_df.to_dict('records')
        main_collection = self.db['sg_geneset_gsea_exp']

        try:
            self.create_db_table('sg_geneset_gsea_exp', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=(gsea_exp), code="53703610")
        else:
            self.bind_object.logger.info("导入gsea_exp：%s出错!" % (gsea_exp))
