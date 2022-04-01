# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.ref_rna_v2.api_base import ApiBase
import pandas as pd
from bson.objectid import ObjectId

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class ExpressPrCorr(ApiBase):
    def __init__(self, bind_object):
        super(ExpressPrCorr, self).__init__(bind_object)
        self._project_type = 'dia'

    @report_check
    def add_expcorr(self, upload_dir, work_dir, main_id=None, corr_way=None,
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
                desc='ExpressPrCorr main table',
                params=params,
                status="start"
            )
            main_id = self.create_db_table('sg_express_pr_corr', [main_info])
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        df_exp = pd.read_table(work_dir + "/express_correlation_info.xls", header=0, sep="\t")
        # if len(df_exp.index) == 0:
        if df_exp.empty:
            self.bind_object.set_error("结果为空，请重新选择参数运行")
        else:
            df_exp = df_exp.round(6)
            df_exp['expcorr_id'] = main_id
            df_exp_dict = df_exp.to_dict("records")
            self.create_db_table('sg_express_pr_corr_detail', df_exp_dict)
            self.update_db_record('sg_express_pr_corr', main_id, status="end",
                                  main_id=main_id, json_dir=upload_dir + "/record.json")

if __name__ == '__main__':
    anno = ExpCorrsf(None)
    work_dir = '/mnt/ilustre/users/sanger-dev/sg-users/litangjian'
    anno.add_ExpCorrsf(work_dir=work_dir, corr_way='spearman', padjust_way='fdr_bh')
