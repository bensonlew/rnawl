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
from mbio.api.database.denovo_rna_v2.api_base import ApiBase

class Gsea(ApiBase):
    def __init__(self, bind_object):
        super(Gsea, self).__init__(bind_object)

    def run_webroot(self, main_id, gsea_dir, term=None):
        task_id = self.bind_object.sheet.id
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或者其对应的字符串！", code="52002401")

        self.bind_object.logger.info("开始导入gsea结果")
        for path in ['gsea_report.xls', 'all_exp.detail', 'all_sets.detail']:
            if os.path.exists(os.path.join(gsea_dir, path)):
                pass
            else:
                self.bind_object.set_error("%s 不存在", variables=(os.path.join(gsea_dir, path)), code="52002402")
        gsea_id = ObjectId(main_id)
        self.add_gsea_stat(gsea_id, os.path.join(gsea_dir, 'gsea_report.xls'))
        self.add_gsea_geneset(gsea_id, os.path.join(gsea_dir, 'all_sets.detail'))
        self.add_gsea_exp(gsea_id, os.path.join(gsea_dir, 'all_exp.detail'))
        main_collection = self.db['sg_geneset_gsea']
        main_collection.update({"_id": main_id},
                               {"$set": {"gsea_sets": term , "status": "end"}})


    @report_check
    def add_gsea_stat(self, gsea_id, gsea_stat):
        gsea_stat_df = pd.read_table(gsea_stat, header=0)
        gsea_stat_df = gsea_stat_df.fillna("")
        gsea_stat_df["gsea_id"] = gsea_id
        row_dict_list = gsea_stat_df.to_dict('records')
        main_collection = self.db['sg_geneset_gsea_stat']

        try:
            self.create_db_table('sg_geneset_gsea_stat', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=(gsea_stat), code="52002403")
        else:
            self.bind_object.logger.info("导入gsea_stat：%s出错!" % (gsea_stat))

    @report_check
    def add_gsea_geneset(self, gsea_id, gsea_geneset):
        gsea_geneset_df = pd.read_table(gsea_geneset, header=0)
        gsea_geneset_df['gsea_id'] = gsea_id
        row_dict_list = gsea_geneset_df.to_dict('records')
        main_collection = self.db['sg_geneset_gsea_geneset']

        try:
            self.create_db_table('sg_geneset_gsea_geneset', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=(gsea_geneset), code="52002404")
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
            self.bind_object.set_error("导入main: %s出错!" , variables=(gsea_exp), code="52002405")
        else:
            self.bind_object.logger.info("导入gsea_exp：%s出错!" % (gsea_exp))
