# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.protein_transcript.api_base import ApiBase
import json
from bson.objectid import ObjectId
import types
import pandas as pd
from scipy.stats import pearsonr
import math


# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class P2gDiff(ApiBase):
    def __init__(self, bind_object):
        super(P2gDiff, self).__init__(bind_object)
        self._project_type = 'itraq_tmt'

    #@report_check
    def add_diff_venn(self, main_id = None, venn_file=None):
        if len(open(venn_file,'rU').readlines())>1:
            venn_data = dict()
            with open(venn_file, 'r') as r_f:
                for line in r_f.readlines():
                    line = line.strip().split('\t')
                    if line:
                        try:
                            venn_data[line[0]] = line[1]
                        except:
                            venn_data[line[0]] = ''
            if not isinstance(main_id, ObjectId):
                if isinstance(main_id, types.StringTypes):
                    main_id = ObjectId(main_id)
                else:
                    raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
            venn_data.update(
                dict(p2g_diff_id = main_id)
            )
            self.create_db_table('sg_p2g_diff_venn', [venn_data])

    def add_diff_nine(self, main_id = None, nine_file=None, fc_file = None):
        if len(open(nine_file,'rU').readlines())>1:
            result_pd = pd.read_table(nine_file, index_col=False)
            result_pd = result_pd.fillna('_')
            corr_value, corr_p = pearsonr(result_pd['protein_fc'], result_pd['transcript_fc'])
            result_pd['p2g_diff_id'] = ObjectId(main_id)
            result_dict_list = result_pd.to_dict("records")
            self.create_db_table("sg_p2g_diff_nine", result_dict_list)
            update_dict = {'corr_value': corr_value, 'corr_p': corr_p, 'main_id': ObjectId(main_id)}
            with open(fc_file) as fc_r:
                for line in fc_r:
                    line = line.strip().split('\t')
                    if len(line) > 0:
                        update_dict.update({line[0]:math.log(float(line[1]), 2)})
            self.update_db_record('sg_p2g_diff', main_id, **update_dict)


