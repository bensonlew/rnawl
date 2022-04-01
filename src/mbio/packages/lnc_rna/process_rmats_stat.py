# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang, shicaiping, qinjincheng'

import re
import os
import subprocess
from collections import defaultdict
import pandas
import numpy as np
import time
from optparse import OptionParser

parser = OptionParser(description='Process output directory of rMATS')
parser.add_option('-i', '--root', dest='root', help='Input result directory')
parser.add_option('-m', '--method', dest='method', choices=['fdr', 'pvalue'], help='Significance test method')
parser.add_option('-c', '--cutoff', dest='cutoff', type=float, help='Hypothesis test cutoff')
parser.add_option('-p', '--psi', dest='psi', help=' Lower limit of Abs of percent spliced in')
(opts, args) = parser.parse_args()

global LEGAL_EVENT_TYPE
global EVENT_FILES
global MATS_FILES

LEGAL_EVENT_TYPE = {'A3SS', 'A5SS', 'MXE', 'RI', 'SE'}
OUTPUT_FILES = {
    'fromGTF.A3SS.txt',
    'fromGTF.A5SS.txt',
    'fromGTF.MXE.txt',
    'fromGTF.novelEvents.A3SS.txt',
    'fromGTF.novelEvents.A5SS.txt',
    'fromGTF.novelEvents.MXE.txt',
    'fromGTF.novelEvents.RI.txt',
    'fromGTF.novelEvents.SE.txt',
    'fromGTF.RI.txt',
    'fromGTF.SE.txt',
    'A3SS.MATS.JC.txt',
    'A3SS.MATS.JCEC.txt',
    'A5SS.MATS.JC.txt',
    'A5SS.MATS.JCEC.txt',
    'MXE.MATS.JC.txt',
    'MXE.MATS.JCEC.txt',
    'RI.MATS.JC.txt',
    'RI.MATS.JCEC.txt',
    'SE.MATS.JC.txt',
    'SE.MATS.JCEC.txt'
}
EVENT_DESC_ITEMS = {
    'A3SS': ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    'A5SS': ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    'SE': ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    'MXE': ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    'RI': ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]
}

def main(root, method, cutoff, psi):
    print 'INFO: start processing {}'.format(root)
    process_rmats_stat(
        root=root,
        pvalue_fdr=method,
        fdr=float(cutoff),
        psi=float(psi)
    )
    print 'INFO: succeed in processing {}'.format(root)

def process_rmats_stat(root, pvalue_fdr, fdr, psi):
    stat_psi = os.path.join(root, 'psi_stats.file.txt')
    events_stats_file = os.path.join(root, 'event_stats.file.txt')

    print('开始获取事件基本信息统计')
    tmp_files = os.listdir(root)
    files = []
    for file in tmp_files:
        files.append(os.path.join(root, file))
    event_info_dic = get_event_stats(files=files, pvalue_fdr=pvalue_fdr, fdr=fdr, psi=psi)

    print('开始写事件统计文件')
    write_event_stat_file(event_stat_f=events_stats_file, event_stat_dic=event_info_dic)

    print('开始读psi info dic')
    psi_info_dic = get_event_psi_dic(event_info_dic)
    write_psi_detail_file(file_src_dic=event_info_dic)
    write_psi_stats_file(stat_file=stat_psi, psi_dic=psi_info_dic)

def get_event_stats(files, pvalue_fdr="fdr", fdr=0.05, psi=0):
    d = dict.fromkeys(LEGAL_EVENT_TYPE)
    for k in LEGAL_EVENT_TYPE:
        d[k] = {
                'JunctionCountOnly_event_id_set': set(), 'ReadsOnTargetAndJunctionCounts_event_id_set': set(),
                'JunctionCountOnly_event_id_set_no': 0, 'ReadsOnTargetAndJunctionCounts_event_id_set_no': 0,

                'JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set': set(),
                'JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no': 0,
                'JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set': set(),
                'JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no': 0,
                }
    for f in files:
        mats_m = re.match(r'(\S+?)\.MATS.(JC|JCEC)\.alter_id\.txt', os.path.basename(f))
        if mats_m:
            event_type = mats_m.group(1)
            print('要读的文件是: %s' % f)
            d[event_type][mats_m.group(2) + '_file'] = f
            data = pandas.read_table(f, sep='\t')
            if pvalue_fdr.lower() == "fdr":
                data_filter1 = data[data['FDR']<=fdr]
                data_filter = data_filter1[abs(data_filter1['IncLevelDifference']) >= psi]
            else:
                data_filter1 = data[data['PValue']<=fdr]
                data_filter = data_filter1[abs(data_filter1['IncLevelDifference']) >= psi]
            if mats_m.group(2) == "JC":
                d[event_type]['JunctionCountOnly_event_id_set_no'] = len(set(data_filter['ID']))
                d[event_type]['JunctionCountOnly_event_id_set'] = set(data_filter['ID'])
            elif mats_m.group(2) == "JCEC":
                d[event_type]['ReadsOnTargetAndJunctionCounts_event_id_set_no'] = len(set(data_filter['ID']))
                d[event_type]['ReadsOnTargetAndJunctionCounts_event_id_set'] = set(data_filter['ID'])
            else:
                d[event_type][mats_m.group(2) + '_event_id_set_no'] = len(set(data_filter['ID']))
                d[event_type][mats_m.group(2) + '_event_id_set'] = set(data_filter['ID'])
            continue

    for as_type in LEGAL_EVENT_TYPE:
        d[as_type]['JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set'] = d[as_type][
                                                                                     'JunctionCountOnly_event_id_set'] & \
                                                                                 d[as_type][
                                                                                     'ReadsOnTargetAndJunctionCounts_event_id_set']
        d[as_type]['JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no'] = len(
            d[as_type]['JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set'])
        d[as_type]['JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set'] = d[as_type][
                                                                                    'JunctionCountOnly_event_id_set'] | \
                                                                                d[as_type][
                                                                                    'ReadsOnTargetAndJunctionCounts_event_id_set']
        d[as_type]['JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no'] = len(
            d[as_type]['JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set'])

    d['total_JunctionCountOnly_event_id_set_no'] = sum(
        [d[e]['JunctionCountOnly_event_id_set_no'] for e in LEGAL_EVENT_TYPE])
    d['total_JunctionCountOnly_event_id_set'] = union_set(
        [d[e]['JunctionCountOnly_event_id_set'] for e in LEGAL_EVENT_TYPE])
    d['total_ReadsOnTargetAndJunctionCounts_event_id_set_no'] = sum(
        [d[e]['ReadsOnTargetAndJunctionCounts_event_id_set_no'] for e in LEGAL_EVENT_TYPE])
    d['total_ReadsOnTargetAndJunctionCounts_event_id_set'] = union_set(
        [d[e]['ReadsOnTargetAndJunctionCounts_event_id_set'] for e in LEGAL_EVENT_TYPE])

    return d

def write_event_stat_file(event_stat_dic, event_stat_f):
    fw = open(event_stat_f, 'wb')
    fw.write('stat_item\tvalue\n')
    for info_key in event_stat_dic.keys():
        if info_key.strip() in LEGAL_EVENT_TYPE:
            for type_info_key in event_stat_dic[info_key].keys():
                if re.match(r'.+_no$', type_info_key.strip()):
                    fw.write('%s\t%d\n' % (info_key + "_" + type_info_key, event_stat_dic[info_key][type_info_key]))
                    continue
        if re.match(r'.+_no$', info_key.strip()):
            fw.write('%s\t%d\n' % (info_key, event_stat_dic[info_key]))
            continue
    fw.close()

def get_event_psi_dic(d, pvalue_fdr="fdr", fdr=0.05, psi=0):
    psi_dic = dict.fromkeys(LEGAL_EVENT_TYPE)
    psi_dic['SAMPLE_1'] = {'JunctionCountOnly': {'exclusion_total': 0, 'inclusion_total': 0, 'total': 0},
                           'ReadsOnTargetAndJunctionCounts': {'exclusion_total': 0, 'inclusion_total': 0, 'total': 0}}
    psi_dic['SAMPLE_2'] = {'JunctionCountOnly': {'exclusion_total': 0, 'inclusion_total': 0, 'total': 0},
                           'ReadsOnTargetAndJunctionCounts': {'exclusion_total': 0, 'inclusion_total': 0, 'total': 0}}
    psi_dic['total'] = {'JunctionCountOnly': 0, 'ReadsOnTargetAndJunctionCounts': 0}
    for k in LEGAL_EVENT_TYPE:
        psi_dic[k] = {'JunctionCountOnly_data': {},
                      # data字典里格式为{e_id1:(pvalue,fdr,inclevel_1,inclevel_2,sample1-sample_level_2,sample_level_2-sample_level_2)}
                      'ReadsOnTargetAndJunctionCounts_data': {},
                      'SAMPLE_1': {'JunctionCountOnly': {'exclusion': 0, 'inclusion': 0, 'total': 0},
                                   'ReadsOnTargetAndJunctionCounts': {'exclusion': 0, 'inclusion': 0, 'total': 0}},
                      'SAMPLE_2': {'JunctionCountOnly': {'exclusion': 0, 'inclusion': 0, 'total': 0},
                                   'ReadsOnTargetAndJunctionCounts': {'exclusion': 0, 'inclusion': 0, 'total': 0}}}
    for event_type in LEGAL_EVENT_TYPE:
        print event_type
        print d[event_type]['JC_file']
        print d[event_type]['JCEC_file']
        JC_file = ''
        JCEC_file = ''
        try:
            JC_file = d[event_type]['JC_file']
            JCEC_file = d[event_type]['JCEC_file']
        except Exception as e:
            raise Exception('%s:不全面的信息字典：%s' % (e, d[event_type]))
        finally:
            psi_dic['JC_file'] = JC_file
            psi_dic['JCEC_file'] = JCEC_file

            jc_data = pandas.read_table(JC_file, index_col=0, sep='\t',
                                        dtype={'Pvalue': np.float64, 'FDR': np.float64,
                                               'IncLevelDifference': np.float64})
            tmp_jc = jc_data[['PValue', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference']]
            if pvalue_fdr.lower() == "fdr":
                data_filter1 = tmp_jc[tmp_jc['FDR']<=fdr]
                jc = data_filter1[abs(data_filter1['IncLevelDifference']) >= psi]
            else:
                data_filter1 = tmp_jc[tmp_jc['PValue']<=fdr]
                jc = data_filter1[abs(data_filter1['IncLevelDifference']) >= psi]
            jc_it = jc.iterrows()
            try:
                while True:
                    val = jc_it.next()
                    event_id = val[0]
                    inc_1 = [float(e) for e in str(val[1]['IncLevel1']).split(',') if isnumber(e)]
                    inc_2 = [float(e) for e in str(val[1]['IncLevel2']).split(',') if isnumber(e)]
                    one_minus_two = sum(inc_1) / float(len(inc_1)) - sum(inc_2) / float(len(inc_2))
                    two_minus_one = -one_minus_two
                    psi_dic[event_type]['JunctionCountOnly_data'][event_id] = (
                        val[1]['PValue'], val[1]['FDR'], val[1]['IncLevel1'], val[1]['IncLevel2'], inc_1, inc_2,
                        one_minus_two, two_minus_one)
                    if one_minus_two == 0:
                        continue
                    psi_dic['total']['JunctionCountOnly'] += 1
                    psi_dic['SAMPLE_1']['JunctionCountOnly']['total'] += 1
                    psi_dic['SAMPLE_2']['JunctionCountOnly']['total'] += 1
                    psi_dic[event_type]['SAMPLE_1']['JunctionCountOnly']['total'] += 1
                    psi_dic[event_type]['SAMPLE_2']['JunctionCountOnly']['total'] += 1
                    if one_minus_two > 0:
                        psi_dic['SAMPLE_1']['JunctionCountOnly']['inclusion_total'] += 1
                        psi_dic['SAMPLE_2']['JunctionCountOnly']['exclusion_total'] += 1
                        psi_dic[event_type]['SAMPLE_1']['JunctionCountOnly']['inclusion'] += 1
                        psi_dic[event_type]['SAMPLE_2']['JunctionCountOnly']['exclusion'] += 1
                    else:
                        psi_dic['SAMPLE_2']['JunctionCountOnly']['inclusion_total'] += 1
                        psi_dic['SAMPLE_1']['JunctionCountOnly']['exclusion_total'] += 1
                        psi_dic[event_type]['SAMPLE_2']['JunctionCountOnly']['inclusion'] += 1
                        psi_dic[event_type]['SAMPLE_1']['JunctionCountOnly']['exclusion'] += 1
            except StopIteration:
                pass
            all_data = pandas.read_table(JCEC_file, index_col=0,
                                         sep='\t',
                                         dtype={'PValue': np.float64, 'FDR': np.float64,
                                                'IncLevelDifference': np.float64})
            tmp_all = all_data[['PValue', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference']]
            all = tmp_all[tmp_all['FDR']<=0.05]
            all_it = all.iterrows()
            #all_it = all_data[['PValue', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference']].iterrows()
            try:
                while True:
                    val = all_it.next()
                    event_id = val[0]
                    inc_1 = [float(e) for e in str(val[1]['IncLevel1']).split(',') if isnumber(e)]
                    inc_2 = [float(e) for e in str(val[1]['IncLevel2']).split(',') if isnumber(e)]
                    one_minus_two = sum(inc_1) / float(len(inc_1)) - sum(inc_2) / float(len(inc_2))
                    two_minus_one = -one_minus_two
                    psi_dic[event_type]['ReadsOnTargetAndJunctionCounts_data'][event_id] = (
                        val[1]['PValue'], val[1]['FDR'], val[1]['IncLevel1'], val[1]['IncLevel2'], inc_1, inc_2,
                        one_minus_two, two_minus_one)
                    if one_minus_two == 0:
                        continue
                    psi_dic['total']['ReadsOnTargetAndJunctionCounts'] += 1
                    psi_dic['SAMPLE_1']['ReadsOnTargetAndJunctionCounts']['total'] += 1
                    psi_dic['SAMPLE_2']['ReadsOnTargetAndJunctionCounts']['total'] += 1
                    psi_dic[event_type]['SAMPLE_1']['ReadsOnTargetAndJunctionCounts']['total'] += 1
                    psi_dic[event_type]['SAMPLE_2']['ReadsOnTargetAndJunctionCounts']['total'] += 1
                    if one_minus_two > 0:
                        psi_dic['SAMPLE_1']['ReadsOnTargetAndJunctionCounts']['inclusion_total'] += 1
                        psi_dic['SAMPLE_2']['ReadsOnTargetAndJunctionCounts']['exclusion_total'] += 1
                        psi_dic[event_type]['SAMPLE_1']['ReadsOnTargetAndJunctionCounts']['inclusion'] += 1
                        psi_dic[event_type]['SAMPLE_2']['ReadsOnTargetAndJunctionCounts']['exclusion'] += 1
                    else:
                        psi_dic['SAMPLE_2']['ReadsOnTargetAndJunctionCounts']['inclusion_total'] += 1
                        psi_dic['SAMPLE_1']['ReadsOnTargetAndJunctionCounts']['exclusion_total'] += 1
                        psi_dic[event_type]['SAMPLE_2']['ReadsOnTargetAndJunctionCounts']['inclusion'] += 1
                        psi_dic[event_type]['SAMPLE_1']['ReadsOnTargetAndJunctionCounts']['exclusion'] += 1
            except StopIteration:
                pass
    return psi_dic

def write_psi_detail_file(file_src_dic):
    for as_type in LEGAL_EVENT_TYPE:
        JunctionCountOnly_file = file_src_dic[as_type]['JC_file']
        ReadsOnTargetAndJunctionCounts_file = file_src_dic[as_type]['JCEC_file']
        new_JunctionCountOnly_file = os.path.join(os.path.dirname(JunctionCountOnly_file),
                                                  re.sub(r'(\S+)\.txt', '\g<1>.psi_info.txt',
                                                         os.path.basename(JunctionCountOnly_file)))
        new_ReadsOnTargetAndJunctionCounts_file = os.path.join(os.path.dirname(ReadsOnTargetAndJunctionCounts_file),
                                                               re.sub(r'(\S+)\.txt', '\g<1>.psi_info.txt',
                                                                      os.path.basename(
                                                                          ReadsOnTargetAndJunctionCounts_file)))
        add_psi_info(JunctionCountOnly_file, new_file=new_JunctionCountOnly_file)
        add_psi_info(ReadsOnTargetAndJunctionCounts_file, new_file=new_ReadsOnTargetAndJunctionCounts_file)
        file_src_dic[as_type]['JunctionCountOnly_add_psi_file'] = new_JunctionCountOnly_file
        file_src_dic[as_type]['ReadsOnTargetAndJunctionCounts_add_psi_file'] = new_ReadsOnTargetAndJunctionCounts_file
    return file_src_dic

def write_psi_stats_file(stat_file, psi_dic):
    fw = open(stat_file, 'wb')
    fw.write('stat_item\tvalue\n')
    for info_key in psi_dic.keys():
        if info_key.strip() in LEGAL_EVENT_TYPE:
            for sample in psi_dic[info_key].keys():
                if re.match(r'SAMPLE_\d$', sample.strip()):
                    for data_src in psi_dic[info_key][sample].keys():
                        for value_item in psi_dic[info_key][sample][data_src].keys():
                            fw.write(
                                '%s\t%d\n' % ("_".join([info_key, sample, data_src, value_item]),
                                              psi_dic[info_key][sample][data_src][value_item]))

        if re.match(r'SAMPLE_\d$', info_key.strip()):
            for data_src in psi_dic[info_key].keys():
                for value_item in psi_dic[info_key][data_src].keys():
                    fw.write(
                        '%s\t%d\n' % ("_".join([info_key, data_src, value_item]),
                                      psi_dic[info_key][data_src][value_item]))
        if re.match(r'total$', info_key.strip()):
            for data_src in psi_dic[info_key].keys():
                fw.write(
                    '%s\t%d\n' % ("_".join([info_key, data_src]),
                                  psi_dic[info_key][data_src]))

    fw.close()

def union_set(set_lst):
    uni = set()
    for s in set_lst:
        if isinstance(s, set):
            uni = uni | s
        else:
            raise Exception('这个元素不是集合类型')
    return uni

def isnumber(aString):
    try:
        float(aString)
        return True
    except:
        return False

def add_psi_info(mats_file, new_file):
    data = pandas.read_table(mats_file, sep='\t',
                             dtype={'PValue': float, 'FDR': float, 'IncLevelDifference': float, 'IncFormLen': int,
                                    'SkipFormLen': int, 'ID.1': str})
    fw = open(new_file, 'wb')
    fw.write('{}\t{}\n'.format('\t'.join(data.keys()[0:len(data.keys()) - 1]), '\t'.join(
        ['average_IncLevel1', 'average_IncLevel2', 'IncLevelDifference', 'increase_inclusion_SAMPLE1',
         'increase_exclusion_SAMPLE1',
         'increase_inclusion_SAMPLE2', 'increase_exclusion_SAMPLE2'])))

    it = data.iterrows()
    while 1:
        try:
            record = it.next()
            IncLevel1_arr = [float(e) for e in str(record[1].IncLevel1).split(',') if isnumber(e)]
            IncLevel2_arr = [float(e) for e in str(record[1].IncLevel2).split(',') if isnumber(e)]
            aver_IncLevel1 = sum(IncLevel1_arr) / float(len(IncLevel1_arr))
            aver_IncLevel2 = sum(IncLevel2_arr) / float(len(IncLevel2_arr))
            increase_inclusion_SAMPLE1 = ''
            increase_inclusion_SAMPLE2 = ''
            increase_exclusion_SAMPLE1 = ''
            increase_exclusion_SAMPLE2 = ''
            if float(record[1].IncLevelDifference) < 0:
                increase_inclusion_SAMPLE1 = 'no'
                increase_inclusion_SAMPLE2 = 'yes'
                increase_exclusion_SAMPLE1 = 'yes'
                increase_exclusion_SAMPLE2 = 'no'

            if float(record[1].IncLevelDifference) > 0:
                increase_inclusion_SAMPLE1 = 'yes'
                increase_inclusion_SAMPLE2 = 'no'
                increase_exclusion_SAMPLE1 = 'no'
                increase_exclusion_SAMPLE2 = 'yes'
            if float(record[1].IncLevelDifference) == 0:
                increase_inclusion_SAMPLE1 = 'no_difference'
                increase_inclusion_SAMPLE2 = 'no_difference'
                increase_exclusion_SAMPLE1 = 'no_difference'
                increase_exclusion_SAMPLE2 = 'no_difference'
            newline = '{}\t{}\t{}\t{}\n'.format(
                '\t'.join([str(e) for e in record[1].get_values()[0:len(record[1].get_values()) - 1]]),
                '\t'.join([str(aver_IncLevel1), str(aver_IncLevel2)]),
                str(record[1].get_values()[len(record[1].get_values()) - 1]), '\t'.join(
                    [increase_inclusion_SAMPLE1, increase_exclusion_SAMPLE1, increase_inclusion_SAMPLE2,
                     increase_exclusion_SAMPLE2]))
            fw.write(newline)
        except StopIteration:
            break

    fw.close()

if __name__ == '__main__':
    if opts.root and opts.method and hasattr(opts, 'cutoff') and hasattr(opts, 'psi'):
        main(opts.root, opts.method, opts.cutoff, opts.psi)
    else:
        parser.print_help()