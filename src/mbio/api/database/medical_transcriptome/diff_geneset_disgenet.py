# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
import datetime
import types
import json
import os
from mbio.api.database.medical_transcriptome.api_base import ApiBase
import unittest
import pandas as pd
import math
from biocluster.config import Config
from bson.son import SON
import glob
import re
from collections import OrderedDict



class DiffGenesetDisgenet(ApiBase):
    def __init__(self, bind_object):
        super(DiffGenesetDisgenet, self).__init__(bind_object)
        self._project_type = 'medical_transcriptome'

    @report_check
    def add_annotation_detail(self, enrich_path, g2e_path, disgenet_id,):
        """
        添加DisGeNET富集详情表
        """
        if not isinstance(disgenet_id, ObjectId):
            if isinstance(disgenet_id, types.StringTypes):
                disgenet_id = ObjectId(disgenet_id)
            else:
                self.bind_object.set_error('disgenet_enrich_id必须为ObjectId对象或其对应的字符串!',)

        if not os.path.exists(enrich_path):
            self.bind_object.set_error('disgenet_enrich_result所指定的路径:%s不存在，请检查！', variables=(enrich_path,),)

        disgenet_df = pd.read_table(enrich_path, header=0)
        disgenet_df['seq_list'] = disgenet_df["GeneID"].map(lambda g: g.split("/"))
        disgenet_df['seq_str'] = disgenet_df["GeneID"].map(lambda g: g.replace('/', ";"))

        disgenet_df['enrich_factor'] = disgenet_df["GeneRatio"].map(lambda x: float(x.split("/")[0])) / \
                                       disgenet_df["BgRatio"].map(lambda x: float(x.split("/")[0]))

        pvalues_min = min([p for p in disgenet_df['pvalue'] if p > 0])
        log_pvalues_min = - math.log10(pvalues_min)
        disgenet_df["neg_log10p_uncorrected"] = disgenet_df['pvalue'].map(
            lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

        pvalues_min = min([p for p in disgenet_df['p.adjust'] if p > 0])
        log_pvalues_min = - math.log10(pvalues_min)
        disgenet_df["neg_log10p_corrected"] = disgenet_df['p.adjust'].map(
            lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

        disgenet_df['disgenet_enrich_id'] = disgenet_id
        disgenet_df.rename(columns={'ID': 'disease_id', 'Description': 'disease_name', 'GeneRatio': 'ratio_in_study',
                                    'BgRatio': 'ratio_in_pop', 'p.adjust': 'padjust', 'Count': 'count'}, inplace=True)
        disgenet_df.drop(['qvalue', 'EntrezID', 'GeneName', 'GeneID'], axis=1, inplace=True)
        data_list = disgenet_df.to_dict('records')

        try:
            collection = self.db['sg_diff_geneset_disgenet_enrich_detail']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.info("导入disgenet富集详情表：%s信息出错:%s" % (enrich_path, e))
        else:
            self.bind_object.logger.info("导入%s disgenet富集详情表信息成功!" % enrich_path)
            record_dict = {"_id": disgenet_id}
            self.update_db_record('sg_diff_geneset_disgenet_enrich',
                                  main_id=disgenet_id, status="end", query_dict=record_dict)
            self.bind_object.logger.info("更新%s disgenet富集主表信息成功!" % disgenet_id)

    @report_check
    def update_main_table_status(self, disgenet_id):
        """
        update DisGeNET富集主表 status
        """
        if not isinstance(disgenet_id, ObjectId):
            if isinstance(disgenet_id, types.StringTypes):
                disgenet_id = ObjectId(disgenet_id)
            else:
                self.bind_object.set_error('disgenet_enrich_id必须为ObjectId对象或其对应的字符串!', )
        record_dict = {"_id": disgenet_id}
        self.update_db_record('sg_diff_geneset_disgenet_enrich',
                              main_id=disgenet_id, status="end", query_dict=record_dict)
        self.bind_object.logger.info("更新%s disgenet富集主表信息成功!" % disgenet_id)


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = dict(id='disgenet_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                    type='workflow',
                    name='medical_transcriptome.medical_transcriptome_test_api',
                    options={

            },
        )
        wsheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wsheet)
        wf.sheet.id = 'tsg_medical_transcriptome'
        wf.sheet.project_sn = '188_5d01dede4f911'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("medical_transcriptome.diff_geneset_disgenet")
        wf.test_api.add_main_table(
            enrich_path="/mnt/ilustre/users/sanger-dev/workspace/20200908/GenesetDisgenet_geneset_disgenet_3840_1697/DisgenetEnrich/output/enrichment.txt",
            geneset_id="5f45c70417b2bf78d9c9c174",
            padjust_method="BH",
        )


if __name__ == '__main__':
    unittest.main()


