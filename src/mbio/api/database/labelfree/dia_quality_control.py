#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/23 12:49
@file    : dia_quality_control.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import sys
import os
import pandas as pd
from collections import OrderedDict
import argparse
import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.labelfree.api_base import ApiBase
import pandas as pd
from bson.objectid import ObjectId
from collections import OrderedDict

class DiaQualityControl(ApiBase):
    def __init__(self, bind_object):
        super(DiaQualityControl, self).__init__(bind_object)
        self._project_type = 'labelfree'
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = '_'.join(self.bind_object.sheet.id.split("_")[0:2])

    def add_peplen(self, report_file, params=None,):
        try:
            report_df = pd.read_csv(report_file, sep='\t')
        except:
            report_df = pd.read_excel(report_file, sep='\t')
        peps = list()
        len2num = dict()
        for pep in report_df['PEP.GroupingKey']:
            pep = ''.join(list(filter(str.isalpha, pep)))
            len_ = len(pep)
            if pep in peps or not len_:
                continue
            peps.append(pep)
            if len_ not in len2num:
                len2num[len_] = 0
            len2num[len_] += 1
        del report_df
        self.bind_object.logger.info(len2num)
        len2num_s = OrderedDict()
        for l in sorted(len2num.keys()):
            len2num_s[str(l)] = len2num[l]
        if params is None:
            params = {"software": "DIA"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': self.project_sn,
            'task_id': self.task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'peplen结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'len': len2num_s
            }
        len_id = self.create_db_table('sg_peptide_len', [insert_data])
        self.update_db_record('sg_peptide_len', len_id, status="end",
                              main_id=len_id)
        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            len_file = os.path.join(qc_dir, 'Peptite_length_distribution.xls')
            with open(len_file, 'w') as ew:
                ew.write('length' + '\t' + 'number' + '\n')
                for l, n in len2num_s.items():
                    ew.write(l + '\t' + str(n) + '\n')

    def add_pepnum(self, report_file, params=None,):
        try:
            report_df = pd.read_csv(report_file, sep='\t')
        except:
            report_df = pd.read_excel(report_file, sep='\t')
        pro2peps = dict()
        for pro, pep in zip(report_df['PG.ProteinGroups'], report_df['PEP.GroupingKey']):
            pros = pro.split(';')
            pep = ''.join(list(filter(str.isalpha, pep)))
            for p in pros:
                if not p in pro2peps:
                    pro2peps[p] = list()
                if pep not in pro2peps[p]:
                    pro2peps[p].append(pep)
        del report_df
        num2num = dict()
        for pro, peps in pro2peps.items():
            num = len(peps)
            if num not in num2num:
                num2num[num] = 0
            num2num[num] += 1
        num2num_s = OrderedDict()
        for num in sorted(num2num.keys()):
            num2num_s[str(num)] = num2num[num]
        if params is None:
            params = {"software": "DIA"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': self.project_sn,
            'task_id': self.task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'pepnum结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'num': num2num_s
            }
        num_id = self.create_db_table('sg_peptide_num', [insert_data])
        self.update_db_record('sg_peptide_num', num_id, status="end",
                              main_id=num_id)

        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            num_file = os.path.join(qc_dir, 'Peptite_number_distribution.xls')
            with open(num_file, 'w') as ew:
                ew.write('pep_num' + '\t' + 'number' + '\n')
                for l, n in num2num_s.items():
                    ew.write(l + '\t' + str(n) + '\n')

    def add_error(self, report_file, params=None,):
        try:
            report_df = pd.read_csv(report_file, sep='\t')
        except:
            report_df = pd.read_excel(report_file, sep='\t')
        mz = list()
        for a, b in zip(report_df['FG.PrecMz'], report_df['FG.PrecMzCalibrated']):
            if a and b:
                mz.append([b, (b-a)/b*1000000])
        del report_df
        if params is None:
            params = {"software": "DIA"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': self.project_sn,
            'task_id': self.task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'peperror结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'mz': mz,
            # 'delta': delta
            }
        error_id = self.create_db_table('sg_peptide_error', [insert_data])
        self.update_db_record('sg_peptide_error', error_id, status="end", main_id=error_id)
        qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
        if os.path.exists(qc_dir):
            error_file = os.path.join(qc_dir, 'dMass.xls')
            with open(error_file, 'w') as ew:
                ew.write('m/z [Da]'+'\t'+ 'DeltaM [ppm]' + '\n')
                for m,d in mz:
                    ew.write(str(m)+'\t'+str(d)+'\n')

    def add_proteininfo(self, exp_ana_file, params=None):
        info2num = OrderedDict()
        with open(exp_ana_file, 'r') as er:
            for line in er:
                if not line.strip():
                    continue
                if u'Sparse Profiles - ' in line:
                    info, num = line.strip().split('\t')
                    info = info.split('Sparse Profiles - ')[1]
                    num = num.split(' of ')[0]
                    if ',' in num:
                        num = num.replace(',', '')
                    info2num[info] = int(num)
        new_columns = ['total_spectrum', 'identified_spectrum', 'peptide_num', 'protein_num', 'protein_group_num']
        try:
            info_df = pd.DataFrame([info2num])
            tmp = info_df.columns.tolist()
            tmp[-1], tmp[-2] = tmp[-2], tmp[-1]
            info_df = info_df[tmp]
            info_df.rename(columns=dict(zip(info_df.columns.tolist(), new_columns)), inplace=True)
            data_list = info_df.to_dict('records')
        except:
            data_list = [dict(zip(new_columns, [0,0,0,0,0]))]
        try:
            info_id = self.create_db_table('sg_protein_info', data_list)
            if params is None:
                params = {"software": "DIA"}
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            insert_data = {
                'project_sn': self.project_sn,
                'task_id': self.task_id,
                'params': params if params else "",
                'status': 'start',
                'desc': 'proteininfo结果表',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            }
            self.update_db_record_2('sg_protein_info', info_id, insert_data,
                                    main_id=info_id)
            self.update_db_record('sg_protein_info', info_id, status='end')

            qc_dir = os.path.join(self.bind_object.work_dir, 'qc')
            if os.path.exists(qc_dir):
                try:
                    info_df.to_csv(os.path.join(qc_dir, 'Protein_infomation.xls'), index=False, header=True)
                except:
                    pass
        except:
            print('protein_info 文件不对')

    def add_proteinmw(self, protein_file, params=None):
        try:
            protein_df = pd.read_csv(protein_file, sep='\t')
        except:
            protein_df = pd.read_excel(protein_file, sep='\t')
        mylist = protein_df.loc[:, 'MW [kDa]'].tolist()
        del protein_df
        num_1_20 = num_21_40 = num_41_60 = num_61_80 = num_81_100 = \
            num_101_120 = num_121_140 \
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
            elif i > 301:
                num_300_inf += 1

        mw_dict = OrderedDict()
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
            params = {"software": "DIA"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': self.project_sn,
            'task_id': self.task_id,
            'params': params if params else "",
            'status': 'start',
            'desc': 'proteinmw结果表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
            'len': mw_dict
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

    def add_coverage(self, protein_file, params=None):
        try:
            protein_df = pd.read_csv(protein_file, sep='\t')
        except:
            protein_df = pd.read_excel(protein_file, sep='\t')

        mylist = protein_df.loc[:, 'Coverage [%]'].tolist()
        del protein_df
        num_1 = num_1_5 = num_5_10 = num_10_20 = num_20_40 = \
            num_40_60 = num_60_80 = num_80_inf = 0

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

        mw_dict = OrderedDict()
        mw_dict['<1'] = num_1
        mw_dict['1-5'] = num_1_5
        mw_dict['5-10'] = num_5_10
        mw_dict['10-20'] = num_10_20
        mw_dict['20-40'] = num_20_40
        mw_dict['40-60'] = num_40_60
        mw_dict['60-80'] = num_60_80
        mw_dict['>80'] = num_80_inf
        if params is None:
            params = {"software": "DIA"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        insert_data = {
            'project_sn': self.project_sn,
            'task_id': self.task_id,
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
