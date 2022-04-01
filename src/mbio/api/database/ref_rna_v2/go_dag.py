# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:20200303
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

class GoDag(ApiBase):
    def __init__(self, bind_object):
        super(GoDag, self).__init__(bind_object)

    def run_webroot(self, main_id, relation_result, visual):
        task_id = self.bind_object.sheet.id
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或者其对应的字符串！", code="53703606")

        self.bind_object.logger.info("开始导入go_dag结果")
        self.add_go_dag_detail(main_id, relation_result)

        with open(visual, 'r') as f:
            visual_json = json.loads(f.read())

        self.update_db_record('sg_geneset_go_dag', main_id, visual_map=visual_json)


    @report_check
    def add_go_dag_detail(self, main_id, relation_result):
        gsea_stat_df = pd.read_table(relation_result, header=0)
        gsea_stat_df = gsea_stat_df.fillna("")
        gsea_stat_df["go_dag_id"] = main_id
        row_dict_list = gsea_stat_df.to_dict('records')
        #  main_collection = self.db['sg_geneset_go_dag_detail']

        try:
            self.create_db_table('sg_geneset_go_dag_detail', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入go_dag出错!" )
