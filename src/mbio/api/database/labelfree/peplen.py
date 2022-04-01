# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.labelfree.api_base import ApiBase
import pandas as pd
from collections import Counter
from collections import OrderedDict
import json
import re
import os

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Peplen(ApiBase):
    def __init__(self, bind_object):
        super(Peplen, self).__init__(bind_object)
        self._project_type = 'labelfree'

    #@report_check
    def add_peplen(self, params=None, project_sn='labelfree',
                    task_id='labelfree', main_id=None, peplen_exp=None):

        df = pd.read_table(peplen_exp, sep ="\t", header = 0)
        # num = df.loc[:, 'Peptide'].tolist()
        # num = [re.sub('[^a-zA-Z]','',x) for x in num]
        # mylist = [len(x) for x in num]
        # mylist.sort()
        num = df.loc[:, 'Annotated Sequence'].tolist()
        mylist = [len(x.split(".")[1]) if u'.' in x else len(x) for x in num]
        mylist.sort()
        length_dict =OrderedDict()
        for i in mylist:
            length_dict[str(i)] = mylist.count(i)

        if params is None:
            params = {"software": "peaks 8.5"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'peplen结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'len': length_dict
            }
        num_id = self.create_db_table('sg_peptide_len', [insert_data])
        self.update_db_record('sg_peptide_len', num_id, status="end",
                              main_id=num_id)

        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            len_file = os.path.join(qc_dir, 'Peptite_length_distribution.xls')
            with open(len_file, 'w') as ew:
                ew.write('length' + '\t' + 'number' + '\n')
                for l, n in length_dict.items():
                    ew.write(l + '\t' + str(n) + '\n')

if __name__ == '__main__':
    anno = Peplen(None)
    peplen_exp = '/mnt/ilustre/users/sanger-dev/sg-users/litangjian' \
                '/protein_dev/peptide.xls'
    anno.add_peplen(peplen_exp=peplen_exp)

