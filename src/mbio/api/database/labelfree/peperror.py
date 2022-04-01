# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.labelfree.api_base import ApiBase
import pandas as pd
import json
import os

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Peperror(ApiBase):
    def __init__(self, bind_object):
        super(Peperror, self).__init__(bind_object)
        self._project_type = 'labelfree'

    #@report_check
    def add_peperror(self, params=None, project_sn='labelfree',
                    task_id='labelfree', main_id=None, peperror_exp=None,s3_png_path=None):

        df = pd.read_table(peperror_exp, sep ="\t", header = 0)
        mz = df.loc[:, 'm/z [Da]'].tolist()
        mz = [round(x, 6) for x in mz]
        delta = df.loc[:, 'DeltaM [ppm]'].tolist()
        delta = [round(x, 6) for x in delta]
        mz_delta = zip(mz, delta)
        if params is None:
            params = {"software": "peaks 8.5"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'peperror结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'mz': mz_delta,
            # 'delta': delta
            }
        if s3_png_path:
            insert_data.update({'peperr_png': s3_png_path})
        error_id = self.create_db_table('sg_peptide_error', [insert_data])
        self.update_db_record('sg_peptide_error', error_id, status="end", main_id=error_id)
        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            error_file = os.path.join(qc_dir, 'dMass.xls')
            with open(error_file, 'w') as ew:
                ew.write('m/z [Da]'+'\t'+ 'DeltaM [ppm]' + '\n')
                for m,d in mz_delta:
                    ew.write(str(m) + '\t' + str(d) + '\n')

if __name__ == '__main__':
    anno = Peperror(None)
    peperror_exp = '/mnt/ilustre/users/sanger-dev/sg-users/litangjian' \
                '/protein_dev/psm.xls'
    anno.add_peperror(peperror_exp=peperror_exp)
