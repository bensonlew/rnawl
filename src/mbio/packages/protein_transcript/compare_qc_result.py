#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: fengyitong&fuwenyao


from collections import OrderedDict
import pandas as pd
import os

plant_result = OrderedDict()
plant_result['raw'] = '/mnt/ilustre/users/sanger-dev/workspace/20181219/Single_wenyao_demo_Arabidopsis_thaliana/HiseqReadsStat/output/fastq_stat.xls'
plant_result['qc1'] = '/mnt/ilustre/users/sanger-dev/workspace/20181219/Single_Single_wenyao_testqc-plant15-50-27/HiseqReadsStat/output/fastq_stat.xls'
plant_result['qc2'] = '/mnt/ilustre/users/sanger-dev/workspace/20181219/Single_Single_wenyao_testqc-plant15-52-14/HiseqReadsStat/output/fastq_stat.xls'
plant_result['fastp1'] = '/mnt/ilustre/users/sanger-dev/workspace/20181220/Single_Single_wenyao_testfastq-plant15-50-58/HiseqReadsStat/output/fastq_stat.xls'
plant_result['fastp2'] = '/mnt/ilustre/users/sanger-dev/workspace/20181220/Single_Single_wenyao_testfastq-plant15-51-05/HiseqReadsStat/output/fastq_stat.xls'
plant_result['fastp3'] = '/mnt/ilustre/users/sanger-dev/workspace/20181220/Single_Single_wenyao_testfastq-plant15-51-11/HiseqReadsStat/output/fastq_stat.xls'

animal_result = OrderedDict()
animal_result['raw'] = '/mnt/ilustre/users/sanger-dev/workspace/20181219/Single_wenyao_demo_Mouse/HiseqReadsStat/output/fastq_stat.xls'
animal_result['qc1'] = '/mnt/ilustre/users/sanger-dev/workspace/20181219/Single_Single_wenyao_testqc15-37-02/HiseqReadsStat/output/fastq_stat.xls'
animal_result['qc2'] = '/mnt/ilustre/users/sanger-dev/workspace/20181219/Single_Single_wenyao_testqc15-49-44/HiseqReadsStat/output/fastq_stat.xls'
animal_result['fastp1'] = '/mnt/ilustre/users/sanger-dev/workspace/20181220/Single_Single_wenyao_testfastq15-33-16/HiseqReadsStat/output/fastq_stat.xls'
animal_result['fastp2'] = '/mnt/ilustre/users/sanger-dev/workspace/20181220/Single_Single_wenyao_testfastq15-35-19/HiseqReadsStat/output/fastq_stat.xls'
animal_result['fastp3'] = '/mnt/ilustre/users/sanger-dev/workspace/20181220/Single_Single_wenyao_testfastq16-30-25/HiseqReadsStat/output/fastq_stat.xls'

def melt_tables(order_dict, col, out, prexi):
    column_list = list()
    for ana in order_dict:
        tmp_col = pd.read_table(order_dict[ana], index_col=0, header=0)[col]
        tmp_col.name = ana
        tmp_col.index.name = 'samples'
        column_list.append(tmp_col)
    result_table = pd.concat(column_list, axis=1)
    result_table.to_csv(out + '/' + prexi + col + '.xls', sep='\t')

zhibiao = ['Total_Reads','Total_Bases','Total_Reads_with_Ns','N_Reads%','A%','T%','C%','G%','N%','Error%','Q20%',   'Q30%','GC%']

# for tou in zhibiao:
    # melt_tables(plant_result, tou, '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/test_qc_module', 'plant')
    # melt_tables(animal_result, tou, '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/test_qc_module', 'animal')



animal_log = OrderedDict()
animal_log['qc1'] = 'Single_wenyao_testqc15-37-02'
animal_log['qc2'] = 'Single_wenyao_testqc15-49-44'
animal_log['fastp1'] = 'Single_wenyao_testfastq15-33-16'
animal_log['fastp2'] = 'Single_wenyao_testfastq16-30-25'
animal_log['fastp3'] = 'Single_wenyao_testfastq15-35-19'


plant_log = OrderedDict()
plant_log['qc1'] = 'Single_wenyao_testqc-plant15-50-27'
plant_log['qc2'] = 'Single_wenyao_testqc-plant15-52-14'
plant_log['fastp1'] = 'Single_wenyao_testfastq-plant15-50-58'
plant_log['fastp2'] = 'Single_wenyao_testfastq-plant15-51-05'
plant_log['fastp3'] = 'Single_wenyao_testfastq-plant15-51-11'


def summary_time(order_dict, out, prexi):
    time_sum = OrderedDict()
    def calute_time(time):
        try:
            time = time.split(':')
            time = 3600*int(time[0]) + 60*int(time[1]) + int(time[2])
        except:
            raise
        return time
    for ana in order_dict:
        with open(os.path.join('/mnt/ilustre/users/sanger-dev/workspace/20181219',order_dict[ana],'log.txt'), 'r') as log_r:
            log_info = log_r.readlines()
            start_time = log_info[0].split(' ')[1]
            end_time = log_info[-1].split(' ')[1]
            run_time = str(calute_time(end_time)-calute_time(start_time))+'s'
            time_sum[ana] = run_time
        with open(os.path.join(out, prexi+'time_summary.xls'),'w') as time_w:
            time_w.write('samples\trun_time\n')
            for sam in time_sum:
                time_w.write(sam + '\t' + time_sum[sam] + '\n')

# summary_time(animal_log, '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/test_qc_module', 'animal')
summary_time(plant_log, '/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/test_qc_module', 'plant')