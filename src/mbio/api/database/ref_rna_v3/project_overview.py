# -*- coding: utf-8 -*-
# __author__ = 'fwy 20210628'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.ref_rna_v2.api_base import ApiBase
import re
import pandas as pd
from bson.objectid import ObjectId
import types
import numpy as np
import unittest
from collections import OrderedDict
import traceback


class ProjectOverview(ApiBase):
    def __init__(self, bind_object):
        super(ProjectOverview, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'

    def add_project_overview(self, task_id = None,group_dict=None,sample_num = None,exp_level="T"):
        self.bind_object.logger.info("开始project_overview的导表")
        exp_level = exp_level
        time_now = datetime.datetime.now()
        main_info = dict(
            task_id=task_id,
            created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
            name="project_view_main_table",
            desc='project_view stat main table',
            status="start",
        )
        main_id = self.create_db_table('sg_project_overview', [main_info])

        qc_stat_infos = {}
        mapping_infos = {}
        gene_stat_infos = {}
        #首先获取信息
        group_dict = group_dict
        # group_dict = self.bind_object.option('group_table').prop['group_dict']
        # sample_num = self.bind_object.option('group_table').prop['sample_number']
        sample_num = sample_num
        qc_stat_infos["sample_num"] = sample_num

        #计算质控相关
        qc_table = self.db["sg_specimen"]
        t_result = qc_table.find({"task_id": task_id, "about_qc": "after"})
        qc_df = pd.DataFrame(list(t_result))
        total_bases = round(float(qc_df["total_bases"].sum())/1000/1000/1000,2)
        min_bases = round(float(qc_df["total_bases"].min())/1000/1000/1000,2)
        min_q30 = round(qc_df["q30_rate"].min(),2)
        qc_stat_infos["total_bases"] = total_bases
        qc_stat_infos["min_bases"] = min_bases
        qc_stat_infos["min_q30"] = min_q30

        try:
            #计算比对相关
            mapping_table = self.db["sg_specimen_mapping"]
            m_result = mapping_table.find({"task_id": task_id})
            m_df = pd.DataFrame(list(m_result))
            m_df["mapping_ratio"] = m_df["mapping_reads"].apply(
                lambda x: float(x.split("(")[-1].split(")")[0].split("%")[0]))

            max_ratio = max(m_df["mapping_ratio"])
            min_tatio = min(m_df["mapping_ratio"])
            mapping_infos["max_ratio"] = max_ratio
            mapping_infos["min_ratio"] = min_tatio
        except:
            pass

        try:
            #计算基因类型相关

            samples = list()
            for each in group_dict:
                samples += group_dict[each]
            conn = self.db["sg_exp"]
            conn_detail = self.db["sg_exp_detail"]
            target_cols = OrderedDict(seq_id=1, _id=0,is_new = 1)
            for each in samples:
                target_cols[each] = 1

            #基因层面
            g_exp_main_table = conn.find_one({"task_id":task_id,"exp_level":"G"})
            g_exp_main_id = g_exp_main_table["main_id"]
            exp_records = conn_detail.find({"exp_id": g_exp_main_id}, target_cols)
            g_exp_df = pd.DataFrame(list(exp_records))
            express_g_df = g_exp_df[g_exp_df.loc[:,samples].sum(axis=1)>0]
            gene_express_num = express_g_df.shape[0]
            ref_express_num = express_g_df[~express_g_df["is_new"]].shape[0]
            new_express_num = express_g_df[express_g_df["is_new"]].shape[0]
            gene_stat_infos["gene_express_num"] = gene_express_num
            gene_stat_infos["ref_gene_express_num"] = ref_express_num
            gene_stat_infos["new_gene_express_num"] = new_express_num
            # if self.bind_object.option('level').lower() == 'transcript':
            if exp_level.lower() == 'transcript':
                t_exp_main_table = conn.find_one({"task_id": task_id, "exp_level": "T"})
                t_exp_main_id = t_exp_main_table["main_id"]
                exp_records = conn_detail.find({"exp_id": t_exp_main_id}, target_cols)
                t_exp_df = pd.DataFrame(list(exp_records))
                express_t_df = t_exp_df[t_exp_df.loc[:, samples].sum(axis=1) > 0]
                trans_express_num = express_t_df.shape[0]
                ref_express_num = express_t_df[~express_t_df["is_new"]].shape[0]
                new_express_num = express_t_df[express_t_df["is_new"]].shape[0]
                gene_stat_infos["trans_express_num"] = trans_express_num
                gene_stat_infos["ref_trans_express_num"] = ref_express_num
                gene_stat_infos["new_trans_express_num"] = new_express_num

        except:
            pass
        record_dict = {"_id": main_id, "task_id": task_id}
        try:
            self.update_db_record('sg_project_overview', query_dict=record_dict, qc_stat_infos=qc_stat_infos)
        except:
            estr = traceback.format_exc()
            self.bind_object.logger.info(estr)
        try:
            self.update_db_record('sg_project_overview', query_dict=record_dict,mapping_infos=mapping_infos)
        except:
            estr = traceback.format_exc()
            self.bind_object.logger.info(estr)
        try:
            self.update_db_record('sg_project_overview', query_dict=record_dict, gene_stat_infos = gene_stat_infos)
        except:
            estr = traceback.format_exc()
            self.bind_object.logger.info(estr)
        try:
            self.update_db_record('sg_project_overview', query_dict=record_dict, main_id=main_id, status = "end")
        except:
            estr = traceback.format_exc()
            self.bind_object.logger.info(estr)


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
      from mbio.workflows.ref_rna_v2.refrna_test_api import RefrnaTestApiWorkflow
      from biocluster.wsheet import Sheet
      from biocluster.config import Config
      import random

      data = {
        # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
        "id": "rofq_d2gckg7mkabti7sbrs83j7",
        "project_sn": "rofq_d2gckg7mkabti7sbrs83j7",
        "type": "workflow",
        "name": "ref_rna_v2.refrna_test_api",
        "options": {
        },
      }
      wsheet = Sheet(data=data)
      wf = RefrnaTestApiWorkflow(wsheet)
      wf.IMPORT_REPORT_DATA = True
      wf.IMPORT_REPORT_AFTER_END = False
      wf.test_api = wf.api.api("ref_rna_v3.project_overview")
      wf.config.DBVersion = 0
      Config().DBVerison =0
      project_type = "ref_rna_v2"
      db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
      conn = db["sg_specimen_group"]
      a=conn.find_one({"task_id": "rofq_d2gckg7mkabti7sbrs83j7"})
      group_dict ={}
      #
      # group_dict = {u'CM': [u'CM_1', u'CM_2', u'CM_3', u'CM_4'],
      #  u'Con': [u'Con_1', u'Con_2', u'Con_3', u'Con_4'],
      #  u'Vac': [u'Vac_4', u'Vac_1', u'Vac_2', u'Vac_3']}
      for n, key in enumerate(a["category_names"]):
          group_dict[key] = []
          for s in a["specimen_names"][n]:
              group_dict[key].append(s)
      samples = list()
      for each in group_dict:
          samples += group_dict[each]
      sample_num = len(samples)
      wf.test_api.add_project_overview(task_id = "rofq_d2gckg7mkabti7sbrs83j7",group_dict=group_dict,sample_num = sample_num,exp_level="transcript")

if __name__ == '__main__':
    unittest.main()