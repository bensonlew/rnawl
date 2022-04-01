# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.dia.api_base import ApiBase
import pandas as pd
from collections import OrderedDict
import json
import os

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Proteinmw(ApiBase):
    def __init__(self, bind_object):
        super(Proteinmw, self).__init__(bind_object)
        self._project_type = 'dia'

    #@report_check
    def add_proteinmw(self, params=None, project_sn='dia',
                    task_id='dia', main_id=None, proteinmw_exp=None):

        df = pd.read_table(proteinmw_exp, sep ="\t", header = 0)
        mylist = df.loc[:, 'MW [kDa]'].tolist()
        num_1_20 =  num_21_40 = num_41_60 = num_61_80 = num_81_100 = \
            num_101_120 = num_121_140\
        = num_141_160 = num_161_180 = num_181_200 = num_201_220 = \
            num_221_240 = num_241_260 = num_261_280 = num_281_300 = num_300_inf = 0

        for i in mylist:
            if i <= 21:
                num_1_20 += 1
            elif i <= 41 and i > 21:
                num_21_40 += 1
            elif i <= 61 and i > 41:
                num_41_60 += 1
            elif i <= 81 and i > 61:
                num_61_80 += 1
            elif i <= 101 and i > 81:
                num_81_100 += 1
            elif i <= 121 and i > 101:
                num_101_120 += 1
            elif i <= 141 and i > 121:
                num_121_140 += 1
            elif i <= 161 and i > 141:
                num_141_160 += 1
            elif i <= 181 and i > 161:
                num_161_180 += 1
            elif i <= 201 and i > 181:
                num_181_200 += 1
            elif i <= 221 and i > 201:
                num_201_220 += 1
            elif i <= 241 and i > 221:
                num_221_240 += 1
            elif i <= 261 and i > 241:
                num_241_260 += 1
            elif i <= 281 and i > 261:
                num_261_280 += 1
            elif i <= 301 and i > 281:
                num_281_300 += 1
            elif  i > 301:
                num_300_inf += 1


        mw_dict =OrderedDict()
        mw_dict['1-21'] = num_1_20
        mw_dict['21-41'] = num_21_40
        mw_dict['41-61'] = num_41_60
        mw_dict['61-81'] = num_61_80
        mw_dict['81-101'] = num_81_100
        mw_dict['101-121'] = num_101_120
        mw_dict['121-141'] = num_121_140
        mw_dict['141-161'] = num_141_160
        mw_dict['161-181'] = num_161_180
        mw_dict['181-201'] = num_181_200
        mw_dict['201-221'] = num_201_220
        mw_dict['221-241'] = num_221_240
        mw_dict['241-261'] = num_241_260
        mw_dict['261-281'] = num_261_280
        mw_dict['281-301'] = num_281_300
        mw_dict['>301'] = num_300_inf
        if params is None:
            params = {"software": "Proteome Discoverer"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'proteinmw结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'len': mw_dict,
            'version':'v3'
            }
        mw_id = self.create_db_table('sg_protein_mw', [insert_data])
        self.update_db_record('sg_protein_mw', mw_id, status="end",
                              main_id=mw_id)

        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            mw_file = os.path.join(qc_dir, 'Protein_mw_distribution.xls')
            with open(mw_file, 'w') as ew:
                ew.write('moleweight' + '\t' + 'number' + '\n')
                for mw, n in mw_dict.items():
                    ew.write(mw + '\t' + str(n) + '\n')

if __name__ == '__main__':
    anno = Proteinmw(None)
    proteinmw_exp = '/mnt/ilustre/users/sanger-dev/sg-users/litangjian' \
                '/protein_dev/protein.xls'
    anno.add_proteinmw(proteinmw_exp=proteinmw_exp)

