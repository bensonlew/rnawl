# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.itraq_and_tmt.api_base import ApiBase
import pandas as pd
from collections import OrderedDict
import json
import os

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Coverage(ApiBase):
    def __init__(self, bind_object):
        super(Coverage, self).__init__(bind_object)
        self._project_type = 'itraq_tmt'

    #@report_check
    def add_coverage(self, params=None, project_sn='itraq_tmt',
                    task_id='itraq_tmt', coverage_exp=None,main_id=None):

        df = pd.read_table(coverage_exp, sep ="\t", header = 0)
        try:
            mylist = df.loc[:, 'Coverage'].tolist()
        except:
            mylist = df.loc[:, 'Coverage [%]'].tolist()
        num_1 =  num_1_5 = num_5_10 = num_10_20 = num_20_40 = \
            num_40_60 = num_60_80  = num_80_inf = 0

        for i in mylist:
            if i <= 1:
                num_1 += 1
            elif i <= 5 and i > 1:
                num_1_5 += 1
            elif i <= 10 and i > 5:
                num_5_10 += 1
            elif i <= 20 and i > 10:
                num_10_20 += 1
            elif i <= 40 and i > 20:
                num_20_40 += 1
            elif i <= 60 and i > 40:
                num_40_60 += 1
            elif i <= 80 and i > 60:
                num_60_80 += 1
            elif i > 80:
                num_80_inf += 1

        mw_dict =OrderedDict()
        mw_dict['<1'] = num_1
        mw_dict['1-5'] = num_1_5
        mw_dict['5-10'] = num_5_10
        mw_dict['10-20'] = num_10_20
        mw_dict['20-40'] = num_20_40
        mw_dict['40-60'] = num_40_60
        mw_dict['60-80'] = num_60_80
        mw_dict['>80'] = num_80_inf
        if params is None:
            params = {"software": "Proteome Discoverer"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'coverage结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'len': mw_dict
            }
        mw_id = self.create_db_table('sg_protein_coverage', [insert_data])
        self.update_db_record('sg_protein_coverage', mw_id, status="end",
                              main_id=mw_id)
        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            cover_file = os.path.join(qc_dir, 'Protein_seq_cover_distribution.xls')
            with open(cover_file, 'w') as ew:
                ew.write('coverage' + '\t' + 'number' + '\n')
                for c, n in mw_dict.items():
                    ew.write(c + '\t' + str(n) + '\n')

if __name__ == '__main__':
    anno = Coverage(None)
    coverage_exp = '/mnt/ilustre/users/sanger-dev/sg-users/litangjian' \
                '/protein_dev/protein.xls'
    anno.add_coverage(coverage_exp=coverage_exp)

