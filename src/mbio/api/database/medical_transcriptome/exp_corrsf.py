# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.medical_transcriptome.api_base import ApiBase
import re
import os
import glob
import pandas as pd
from bson.objectid import ObjectId
from collections import OrderedDict
import unittest
from itertools import islice

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class ExpCorrsf(ApiBase):
    def __init__(self, bind_object):
        super(ExpCorrsf, self).__init__(bind_object)
        self._project_type = 'medical_transcriptome'

    @report_check
    def add_ExpCorrsf(self, upload_dir, work_dir, main_id=None, corr_way=None,
                      padjust_way=None, params=None, project_sn=None,
                      task_id=None):
        if main_id is None:
            # prepare main table info
            name = "ExpCoranalysis" + '_' + corr_way + '_' + padjust_way + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                json_dir = work_dir + "/record.json",
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='ExpCorrsf main table',
                params=params,
                status="start"
            )
            main_id = self.create_db_table('sg_exp_corrsf', [main_info])
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        df_exp = pd.read_table(work_dir + "/express_correlation_info.xls", header=0, sep="\t")
        df_exp = df_exp.round(6)
        df_exp['expcorr_id'] = main_id
        df_exp_dict = df_exp.to_dict("records")
        self.create_db_table('sg_exp_corrsf_detail', df_exp_dict)
        self.update_db_record('sg_exp_corrsf', main_id, status="end",
                              main_id=main_id, json_dir=upload_dir + "/record.json")



# if __name__ == '__main__':
#     anno = ExpCorrsf(None)
#     work_dir = '/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong'
#     anno.add_ExpCorrsf(work_dir=work_dir, corr_way='spearman', padjust_way='fdr_bh')



class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
      from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
      from biocluster.wsheet import Sheet
      import random
      data = dict(
          id= "medical_transcriptome_corrsf",
          project_sn= "medical_transcriptome_corrsf",
          type="workflow",
          name="medical_transcriptome.medical_transcriptome_test_api",
                  options={
                  },
                  )
      wsheet = Sheet(data=data)
      wf = MedicalTranscriptomeTestApiWorkflow(wsheet)
      wf.IMPORT_REPORT_DATA = True
      wf.IMPORT_REPORT_AFTER_END = False
      wf.test_api = wf.api.api("medical_transcriptome.exp_corrsf")
      wf.test_api.add_ExpCorrsf(work_dir = '/mnt/ilustre/users/sanger-dev/workspace/20200826/ExpCorrsf_exp_corrsf_workflow_9112_7256/ExpCorrsf/output',
                                upload_dir="/mnt/ilustre/users/sanger-dev/workspace/20200826/ExpCorrsf_exp_corrsf_workflow_9112_7256/ExpCorrsf/output",
                                corr_way='spearman', padjust_way="fdr")
          # upload_dir="/mnt/ilustre/users/sanger-dev/workspace/20200826/ExpCorrsf_exp_corrsf_workflow_9112_7256/ExpCorrsf/output",
          # work_dir= "/mnt/ilustre/users/sanger-dev/workspace/20200826/ExpCorrsf_exp_corrsf_workflow_9112_7256/ExpCorrsf/output")


if __name__ == '__main__':
    unittest.main()