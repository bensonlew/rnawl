# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.labelfree.api_base import ApiBase
import re
import os
import pandas as pd
from bson.objectid import ObjectId
from collections import OrderedDict

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Pca(ApiBase):
    def __init__(self, bind_object):
        super(Pca, self).__init__(bind_object)
        self._project_type = 'labelfree'

    #@report_check
    def add_pca(self, exp_output_dir, exp_id=None, params=None,
                    project_sn='labelfree', task_id='labelfree', main_id=None):
        if main_id is None:
            # prepare main table info
            name = "ExpPCA" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='PCA main table',
                params=params,
                status="start",
                version='v3',
            )
            main_id = self.create_db_table('sg_express_pca', [main_info])
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        # prepare detail table info
        # result, pc_ratio_dict = self.pca(all_exp_pd)
        target_file = os.path.join(exp_output_dir, 'pca_sites.xls')
        result = pd.read_table(target_file, header=0)
        result.rename(columns={"Sample_ID": "sample"}, inplace=True)
        result = result.round(6)
        row_dict_list = result.to_dict('records')
        #
          # result, pc_ratio_dict = self.pca(all_exp_pd)
        target_file = os.path.join(exp_output_dir, 'pca_rotation.xls')
        result = pd.read_table(target_file, header=0)
        result.rename(columns={"Sample_ID": "accession_id"}, inplace=True)
        result = result.round(6)
        gene_dict_list = result.to_dict('records')
        #
        target_file = os.path.join(exp_output_dir, 'pca_importance.xls')
        t = pd.read_table(target_file, header=0)
        t = t.round(6)
        pc_ratio_dict = OrderedDict(zip(list(t.iloc[:, 0]), list(t.iloc[:,1])))
        self.update_db_record('sg_express_pca', main_id, ratio_dict=pc_ratio_dict)
        # insert detail
        tag_dict = dict(pca_id=main_id)
        self.create_db_table('sg_express_pca_sample_detail', row_dict_list,
                             tag_dict=tag_dict)
        self.create_db_table('sg_express_pca_gene_detail', gene_dict_list,
                             tag_dict=tag_dict)
        self.update_db_record('sg_express_pca', main_id, status="end",
                              main_id=main_id, pcn=list(t.iloc[:, 0]))

if __name__ == '__main__':
    anno = Pca(None)
    exp_output_dir = '/mnt/ilustre/users/sanger-dev/workspace/20180308/Single_test_pca40852/Pca/output'
    anno.add_pca(exp_output_dir=exp_output_dir)
