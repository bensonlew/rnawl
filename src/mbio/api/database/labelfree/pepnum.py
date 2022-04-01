# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.labelfree.api_base import ApiBase
import pandas as pd
from collections import OrderedDict
import json
import os

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Pepnum(ApiBase):
    def __init__(self, bind_object):
        super(Pepnum, self).__init__(bind_object)
        self._project_type = 'labelfree'

    #@report_check
    def add_pepnum(self, params=None, project_sn='labelfree',
                    task_id='labelfree', main_id=None, pepnum_exp=None):

        df = pd.read_table(pepnum_exp, sep ="\t", header = 0)
        # num = dict(df.loc[:, '#Peptides'].value_counts())
        num = OrderedDict(df.loc[:, '# Peptides'].value_counts(ascending=True).sort_index())
        keys_1 = [str(int(x)) for x in num.keys()]
        values_1 = num.values()
        num = OrderedDict(zip(keys_1, num.values()))
        if params is None:
            params = {"software": "peaks 8.5"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'pepnum结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'num': num
            }
        num_id = self.create_db_table('sg_peptide_num', [insert_data])
        self.update_db_record('sg_peptide_num', num_id, status="end",
                              main_id=num_id)

        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            num_file = os.path.join(qc_dir, 'Peptite_number_distribution.xls')
            with open(num_file, 'w') as ew:
                ew.write('pep_num' + '\t' + 'number' + '\n')
                for l, n in num.items():
                    ew.write(l + '\t' + str(n) + '\n')

if __name__ == '__main__':
    anno = Pepnum(None)
    pepnum_exp = '/mnt/ilustre/users/sanger-dev/sg-users/litangjian' \
                '/protein_dev/protein.xls'
    anno.add_pepnum(pepnum_exp=pepnum_exp)

