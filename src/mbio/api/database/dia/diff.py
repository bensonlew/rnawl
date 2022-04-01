# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.dia.api_base import ApiBase
import re
import os
import glob
import pandas as pd
from bson.objectid import ObjectId
from collections import OrderedDict
import numpy as np

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Diff(ApiBase):
    def __init__(self, bind_object):
        super(Diff, self).__init__(bind_object)
        self._project_type = 'dia'

    #@report_check
    def add_diff(self, work_dir, compare_dict_xls, num_summary_xls,
        allsummary_xls, group_dict=None, group_id=None,project_sn=None, task_id=None,
        main_id=None,params=None, method_type=None, fill_type=None):
        df_cmp = pd.read_table(compare_dict_xls, header=0, sep="\t")
        cmp_list = list(df_cmp.columns)
        df_cmp_dict = df_cmp.to_dict("list")
        for i in df_cmp_dict.values():
            while np.nan in i:
                i.remove(np.nan)
        df_num = pd.read_table(num_summary_xls, header=0, sep="\t")
        df_num.columns = [x.replace("_vs_", "|") for x in list(df_num.columns)]
        df_num_dict = df_num.to_dict("list")
        if main_id is None:
            # prepare main table info
            name = "Diff" + '_' + method_type + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Diff main table',
                params=params,
                cmp_detail=df_cmp_dict,
                sig_status=df_num_dict,
                cmp_list=cmp_list,
                group_dict=group_dict,
                group_id=group_id,
                # diff_method=method_type,
                status="start",
                version='v3'

            )
            # print('双击666666')
            # print(main_info)
            # print('双击666666')
            main_id = self.create_db_table('sg_diff', [main_info])
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:main_id = ObjectId(main_id)
        # prepare detail table info
        # result, pc_ratio_dict = self.diff(all_exp_pd)
        tag_dict = dict(diff_id=main_id)
        files = glob.glob(work_dir + '/diffcmp*.csv')
        #找到workdor下面所有的diffcmp字段的文件，导入详情表数据,路径需要拼接，不能直接(r"work_dir +
        # "/diffcmp*.csv"
        for file in files:
            result = pd.read_table(file, header=0, sep="\t", dtype = {'accession_id' : str})
            result = result.round(6)
            dict_list = result.to_dict('records')
            self.create_db_table('sg_diff_detail', dict_list, tag_dict=tag_dict)

        # 导入yes_no表的相关数据、

        df_summary = pd.read_table(allsummary_xls, header=0, sep="\t", dtype = {'accession_id' : str})
        df_summary = df_summary.fillna('no')
        len_df = df_summary.shape[1] - 1

        if fill_type.lower() != 'all':
            num_cols = 7
        else:
            num_cols = 4    # ncols per comparison in allsummary_xls

        def num_sum(row):
            count_sum_yes = 0
            for i in range(2, len_df, num_cols):
                if row[i] == "yes":
                    count_sum_yes += 1
            return count_sum_yes

        def num_sum_up(row):
            count_sum_up = 0
            for i in range(3, len_df, num_cols):
                if row[i] == "yes":
                    count_sum_up += 1
            return count_sum_up

        def num_sum_down(row):
            count_sum_down = 0
            for i in range(4, len_df+1, num_cols):
                if row[i] == "yes":
                    count_sum_down += 1
            return count_sum_down


        df_summary["sum"] = df_summary.apply(num_sum, axis=1)
        df_summary["sum_up"] = df_summary.apply(num_sum_up, axis=1)
        df_summary["sum_down"] = df_summary.apply(num_sum_down, axis=1)
        dict_summary = df_summary.to_dict('records')
        self.create_db_table('sg_diff_summary', dict_summary, tag_dict=tag_dict)
        self.update_db_record('sg_diff', main_id, status="end",
                              main_id=main_id, cmp_list=cmp_list, cmp_detail=df_cmp_dict,sig_status=df_num_dict,)

if __name__ == '__main__':
    anno = Diff(None)
    work_dir = '/mnt/ilustre/users/sanger-dev/workspace/20180419/Single_diff9475fyt/Diff'
    compare_dict_xls = '/mnt/ilustre/users/sanger-dev/workspace/20180419/Single_diff9475fyt/Diff/compare_dict.xls'
    num_summary_xls = '/mnt/ilustre/users/sanger-dev/workspace/20180419/Single_diff9475fyt/Diff/num_summary.xls'
    allsummary_xls = '/mnt/ilustre/users/sanger-dev/workspace/20180419/Single_diff9475fyt/Diff/all_labelfree_summary.xls'
    anno.add_diff(work_dir=work_dir, compare_dict_xls=compare_dict_xls,
               num_summary_xls=num_summary_xls, allsummary_xls=allsummary_xls )
