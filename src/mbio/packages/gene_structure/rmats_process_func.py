# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/4/19 08:49

import re, os, Bio, argparse, sys, fileinput, urllib2
import subprocess
from collections import defaultdict
import pandas
import numpy as np
from bson.son import SON
import bson.binary
from cStringIO import StringIO
import time

global LEGAL_EVENT_TYPE
global PRIMARY_DIRS
global EVENT_FILES
global MATS_FILES

LEGAL_EVENT_TYPE = {'A3SS', 'A5SS', 'MXE', 'RI', 'SE'}
PRIMARY_DIRS = {'ASEvents', 'commands.txt', 'config.txt', 'MATS_output', 'SAMPLE_1', 'SAMPLE_2'}
EVENT_FILES = {'fromGTF.A3SS.txt', 'fromGTF.A5SS.txt', 'fromGTF.MXE.txt', 'fromGTF.novelEvents.A3SS.txt',
               'fromGTF.novelEvents.A5SS.txt', 'fromGTF.novelEvents.MXE.txt', 'fromGTF.novelEvents.RI.txt',
               'fromGTF.novelEvents.SE.txt', 'fromGTF.RI.txt', 'fromGTF.SE.txt'
               }
MATS_FILES = {'A3SS.MATS.JunctionCountOnly.txt', 'A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt',
              'A5SS.MATS.JunctionCountOnly.txt', 'A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt',
              'MXE.MATS.JunctionCountOnly.txt', 'MXE.MATS.ReadsOnTargetAndJunctionCounts.txt',
              'RI.MATS.JunctionCountOnly.txt', 'RI.MATS.ReadsOnTargetAndJunctionCounts.txt',
              'SE.MATS.JunctionCountOnly.txt', 'SE.MATS.ReadsOnTargetAndJunctionCounts.txt'}
EVENT_DESC_ITEMS = {'A3SS': ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
                    'A5SS': ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
                    'SE': ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
                    'MXE': ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES",
                            "upstreamEE", "downstreamES", "downstreamEE"],
                    'RI': ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES",
                           "downstreamEE"]}


def modify_id_for_txt(infile, outfile):
    f_name_match = re.match(r'^.*?(A3SS|A5SS|MXE|RI|SE).*?\.txt$', os.path.basename(infile).strip())
    if f_name_match:
        tmp = outfile + '.temp'
        content = re.sub(r'\"', '', open(infile).read())
        open(tmp, 'w').write(content)
        event_type = f_name_match.group(1)
        '''awk '{if ($1~/[0-9]+/) printf "A3SS_"$0"\n";  else print $0}' fromGTF.A3SS.txt'''
        # cmd = """awk -F \'\\t\'  \'{if ($1~/[0-9]+/) printf  "%s"$0"\\n"; else print $0}\'   %s > %s""" % (
        #     event_type + "_", infile, outfile)
        cmd = """awk -F \'\\t\'  \'{if ($1~/[0-9]+/) printf  "%s"$0"\\n"; else print $0}\'   %s > %s""" % (
            event_type + "_", tmp, outfile)
        subprocess.call(cmd, shell=True)

    else:
        raise Exception('输入的rmats结果文件名字有误：%s' % (infile))


def check_rmats_out_dir(dirpath):
    if PRIMARY_DIRS.difference(set(os.listdir(dirpath))):
        raise Exception('不完整的一级rmats结果目录结构: %s' % (set(os.listdir(dirpath))))

    if MATS_FILES.difference(set(os.listdir(os.path.join(dirpath, 'MATS_output')))):
        raise Exception('不完整的mats 结果目录结构： %s，  %s' % (
            os.path.join(dirpath, 'MATS_output'), set(os.listdir(os.path.join(dirpath, 'MATS_output')))))
    else:
        # get event id list for filtering purpose. added by gdq
        event_types = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']
        event_id_dict = dict()
        for each in event_types:
            target_file = each + '.MATS.ReadsOnTargetAndJunctionCounts.txt'
            target_file = os.path.join(dirpath, 'MATS_output', target_file)
            with open(target_file) as f:
                f.readline()
                id_list = [int(x.split('\t')[0]) for x in f if x.strip()]
            event_id_dict[each] = id_list

    if EVENT_FILES.difference(set(os.listdir(os.path.join(dirpath, 'ASEvents')))):
        raise Exception('不完整的as events 目录结构： %s，  %s' % (
            os.path.join(dirpath, 'ASEvents'), set(os.listdir(os.path.join(dirpath, 'ASEvents')))))
    else:
        # modify the files in ASEvents
        for each in event_types:
            target_file = 'fromGTF.{}.txt'.format(each)
            target_file2 = 'fromGTF.novelEvents.{}.txt'.format(each)
            target_file = os.path.join(dirpath, 'ASEvents', target_file)
            target_file2 = os.path.join(dirpath, 'ASEvents', target_file2)
            id_list = event_id_dict[each]
            # modify 'fromGTF.{}.txt'
            target_file_table = pandas.read_table(target_file, index_col=0, sep='\s', header=0)
            target_file_table = target_file_table.query("index in @id_list")
            target_file_table.to_csv(target_file, sep='\t')
            # modify fromGTF.novelEvents.{}.txt
            target_file2_table = pandas.read_table(target_file2, index_col=0, sep='\s', header=0)
            target_file2_table = target_file2_table.query("index in @id_list")
            target_file2_table.to_csv(target_file2, sep='\t')

    to_be_altered_id_files = set([os.path.join(os.path.join(dirpath, 'ASEvents'), name) for name in
                                  os.listdir(os.path.join(dirpath, 'ASEvents')) if
                                  re.match(r'^fromGTF\.(novelEvents\.)?(A3SS|A5SS|MXE|RI|SE)\.txt$', name.strip())] +
                                 [os.path.join(os.path.join(dirpath, 'MATS_output'), name) for name in
                                  os.listdir(os.path.join(dirpath, 'MATS_output')) if re.match(
                                     r'^(A3SS|A5SS|MXE|RI|SE)\.MATS\.(ReadsOnTargetAndJunctionCounts|JunctionCountOnly)\.txt$',
                                     name.strip())])
    info_dic = defaultdict(str)
    for f in to_be_altered_id_files:
        alter_id_file = os.path.join(os.path.dirname(f),
                                     re.match(r'^(\S+)\.txt$', os.path.basename(f).strip()).group(1) + '.alter_id.txt')
        # format_line_file = os.path.join(os.path.dirname(f),
        #                            re.match(r'^(\S+)\.txt$', os.path.basename(f).strip()).group(1) + '.format_line.txt')
        # info_dic[f] = (format_line_file, alter_id_file)
        info_dic[f] = alter_id_file
    return info_dic


def isnumber(aString):
    try:
        float(aString)
        return True
    except:
        return False


def union_set(set_lst):
    uni = set()
    for s in set_lst:
        if isinstance(s, set):
            uni = uni | s
        else:
            raise Exception('这个元素不是集合类型')
    return uni


def get_event_stats(files):
    d = dict.fromkeys(LEGAL_EVENT_TYPE)
    for k in LEGAL_EVENT_TYPE:
        d[k] = {'all_event_file': None, 'novel_event_file': None, 'JunctionCountOnly_file': None,
                'ReadsOnTargetAndJunctionCounts_file': None,  # FILE
                'all_event_id_no': 0, 'novel_event_id_no': 0,
                'all_event_id_set': set(), 'novel_event_id_set': set(), 'old_event_id_set': set(),  #

                'JunctionCountOnly_event_id_set': set(), 'ReadsOnTargetAndJunctionCounts_event_id_set': set(),
                'JunctionCountOnly_event_id_set_no': 0, 'ReadsOnTargetAndJunctionCounts_event_id_set_no': 0,

                'JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set': set(),
                'JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no': 0,
                'JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set': set(),
                'JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no': 0,
                }
    for f in files:
        event_m = re.match(r'fromGTF\.(novelEvents\.)?(A3SS|A5SS|MXE|RI|SE)\.alter_id\.txt', os.path.basename(f))
        mats_m = re.match(r'(\S+?)\.MATS.(JunctionCountOnly|ReadsOnTargetAndJunctionCounts)\.alter_id\.txt',
                          os.path.basename(f))
        if event_m:
            event_type = event_m.group(2)
            print('要读的是文件: %s' % f)
            data = pandas.read_table(f, sep='\t')
            if event_m.group(1) is not None:
                d[event_type]['novel_event_file'] = f
                d[event_type]['novel_event_id_no'] = len(set(data['ID']))
                d[event_type]['novel_event_id_set'] = set(data['ID'])
            else:
                d[event_type]['all_event_file'] = f
                d[event_type]['all_event_id_no'] = len(set(data['ID']))
                d[event_type]['all_event_id_set'] = set(data['ID'])
            continue
        if mats_m:
            event_type = mats_m.group(1)
            d[event_type][mats_m.group(2) + '_file'] = f
            data = pandas.read_table(f, sep='\t')
            data_filter = data[data['FDR']<=0.05]
            d[event_type][mats_m.group(2) + '_event_id_set_no'] = len(set(data_filter['ID']))
            d[event_type][mats_m.group(2) + '_event_id_set'] = set(data_filter['ID'])
            continue
        if not (event_m or mats_m):
            print('这个文件： %s 既不是事件文件 也不是 mats 结果文件' % (f))
    for as_type in LEGAL_EVENT_TYPE:
        d[as_type]['old_event_id_set'] = d[as_type]['all_event_id_set'] - d[as_type]['novel_event_id_set']
        d[as_type]['old_event_id_set_no'] = len(d[as_type]['old_event_id_set'])
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

    d['total_as_events_no'] = sum([d[e]['all_event_id_no'] for e in LEGAL_EVENT_TYPE])
    d['total_as_events_id_set'] = union_set([d[e]['all_event_id_set'] for e in LEGAL_EVENT_TYPE])
    d['total_as_novel_events_no'] = sum([d[e]['novel_event_id_no'] for e in LEGAL_EVENT_TYPE])
    d['total_as_novel_events_id_set'] = union_set([d[e]['novel_event_id_set'] for e in LEGAL_EVENT_TYPE])

    d['total_JunctionCountOnly_event_id_set_no'] = sum(
        [d[e]['JunctionCountOnly_event_id_set_no'] for e in LEGAL_EVENT_TYPE])
    d['total_JunctionCountOnly_event_id_set'] = union_set(
        [d[e]['JunctionCountOnly_event_id_set'] for e in LEGAL_EVENT_TYPE])
    d['total_ReadsOnTargetAndJunctionCounts_event_id_set_no'] = sum(
        [d[e]['ReadsOnTargetAndJunctionCounts_event_id_set_no'] for e in LEGAL_EVENT_TYPE])
    d['total_ReadsOnTargetAndJunctionCounts_event_id_set'] = union_set(
        [d[e]['ReadsOnTargetAndJunctionCounts_event_id_set'] for e in LEGAL_EVENT_TYPE])

    return d


def get_event_psi_dic(d):
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
        JunctionCountOnly_file = ''
        ReadsOnTargetAndJunctionCounts_file = ''
        try:
            JunctionCountOnly_file = d[event_type]['JunctionCountOnly_file']
            ReadsOnTargetAndJunctionCounts_file = d[event_type]['ReadsOnTargetAndJunctionCounts_file']
        except Exception as e:
            raise Exception('%s:不全面的信息字典：%s' % (e, d[event_type]))
        finally:
            psi_dic['JunctionCountOnly_file'] = JunctionCountOnly_file
            psi_dic['ReadsOnTargetAndJunctionCounts_file'] = ReadsOnTargetAndJunctionCounts_file

            jc_data = pandas.read_table(JunctionCountOnly_file, index_col=0, sep='\t',
                                        dtype={'Pvalue': np.float64, 'FDR': np.float64,
                                               'IncLevelDifference': np.float64})
            tmp_jc = jc_data[['PValue', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference']]
            jc = tmp_jc[tmp_jc['FDR']<=0.05]
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
            all_data = pandas.read_table(ReadsOnTargetAndJunctionCounts_file, index_col=0,
                                         sep='\t',
                                         dtype={'Pvalue': np.float64, 'FDR': np.float64,
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


def write_psi_detail_file(file_src_dic):
    for as_type in LEGAL_EVENT_TYPE:
        JunctionCountOnly_file = file_src_dic[as_type]['JunctionCountOnly_file']
        ReadsOnTargetAndJunctionCounts_file = file_src_dic[as_type]['ReadsOnTargetAndJunctionCounts_file']
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


def write_big_detail_file(out, big_file):
    fw = open(big_file, 'wb')
    head_lst = ["event_id", "type", "novel", "old", "gene", "gene_symbol", "chr", "strand",
                "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
                "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE",
                "riExonStart_0base", "riExonEnd", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base",
                "2ndExonEnd",
                "diff_JunctionCountOnly", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
                "IncFormLen_JunctionCountOnly", "SkipFormLen_JunctionCountOnly", "PValue_JunctionCountOnly",
                "FDR_JunctionCountOnly", "IncLevel1_JunctionCountOnly", "IncLevel2_JunctionCountOnly",
                "average_IncLevel1_JunctionCountOnly", "average_IncLevel2_JunctionCountOnly",
                "IncLevelDifference_JunctionCountOnly", "increase_inclusion_SAMPLE1_JunctionCountOnly",
                "increase_exclusion_SAMPLE1_JunctionCountOnly", "increase_inclusion_SAMPLE2_JunctionCountOnly",
                "increase_exclusion_SAMPLE2_JunctionCountOnly",
                "diff_ReadsOnTargetAndJunctionCounts", "IC_SAMPLE_1", "SC_SAMPLE_1", "IC_SAMPLE_2", "SC_SAMPLE_2",
                "IncFormLen_ReadsOnTargetAndJunctionCounts", "SkipFormLen_ReadsOnTargetAndJunctionCounts",
                "PValue_ReadsOnTargetAndJunctionCounts", "FDR_ReadsOnTargetAndJunctionCounts",
                "IncLevel1_ReadsOnTargetAndJunctionCounts", "IncLevel2_ReadsOnTargetAndJunctionCounts",
                "average_IncLevel1_ReadsOnTargetAndJunctionCounts", "average_IncLevel2_ReadsOnTargetAndJunctionCounts",
                "IncLevelDifference_ReadsOnTargetAndJunctionCounts",
                "increase_inclusion_SAMPLE1_ReadsOnTargetAndJunctionCounts",
                "increase_exclusion_SAMPLE1_ReadsOnTargetAndJunctionCounts",
                "increase_inclusion_SAMPLE2_ReadsOnTargetAndJunctionCounts",
                "increase_exclusion_SAMPLE2_ReadsOnTargetAndJunctionCounts",
                "diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts",
                "diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts"
                ]
    head = '\t'.join(head_lst) + "\n"
    fw.write(head)
    for type in LEGAL_EVENT_TYPE:
        line_18_items_dic = dict(zip(head_lst[8:26], range(0, 18)))
        line_18_items = ['null'] * 18
        all_event_file_type = os.path.join(out, 'ASEvents/fromGTF.%s.alter_id.txt' % (type))
        novel_event_file_type = os.path.join(out, 'ASEvents/fromGTF.novelEvents.%s.alter_id.txt' % (type))
        JunctionCountOnly_file_type = os.path.join(out,
                                                   'MATS_output/%s.MATS.JunctionCountOnly.alter_id.psi_info.txt' % (
                                                       type))
        ReadsOnTargetAndJunctionCounts_file_type = os.path.join(out,
                                                                'MATS_output/%s.MATS.ReadsOnTargetAndJunctionCounts.alter_id.psi_info.txt' % (
                                                                    type))
        all_event_info_type = pandas.read_table(all_event_file_type, index_col=0, sep='\t')
        all_event_set_type = set(all_event_info_type.index.tolist())
        novel_event_set_type = set(
            pandas.read_table(novel_event_file_type, index_col=0, sep='\t').index.tolist())
        JunctionCountOnly_info_type = pandas.read_table(JunctionCountOnly_file_type, sep='\t',
                                                        dtype={'FDR': str, 'PValue': str, 'IncLevelDifference': str,
                                                               'average_IncLevel2': str, 'average_IncLevel1': str},
                                                        index_col=0)
        JunctionCountOnly_event_set_type = set(JunctionCountOnly_info_type.index.tolist())
        ReadsOnTargetAndJunctionCounts_info_type = pandas.read_table(ReadsOnTargetAndJunctionCounts_file_type,
                                                                     dtype={'FDR': str, 'PValue': str,
                                                                            'IncLevelDifference': str,
                                                                            'average_IncLevel2': str,
                                                                            'average_IncLevel1': str}, sep='\t',
                                                                     index_col=0)
        ReadsOnTargetAndJunctionCounts_event_set_type = set(ReadsOnTargetAndJunctionCounts_info_type.index.tolist())

        for event_id in all_event_set_type:
            novel = 'no'
            old = 'yes'
            diff_ReadsOnTargetAndJunctionCounts = 'no'
            diff_JunctionCountOnly = 'no'
            diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts = 'no'
            diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts = 'no'

            if event_id in novel_event_set_type:
                novel = 'yes'
                old = 'no'
            if event_id in ReadsOnTargetAndJunctionCounts_event_set_type:
                diff_ReadsOnTargetAndJunctionCounts = 'yes'
            if event_id in JunctionCountOnly_event_set_type:
                diff_JunctionCountOnly = 'yes'
            if (event_id in JunctionCountOnly_event_set_type) and (
                        event_id in ReadsOnTargetAndJunctionCounts_event_set_type):
                diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts = 'yes'
            if (event_id in JunctionCountOnly_event_set_type) or (
                        event_id in ReadsOnTargetAndJunctionCounts_event_set_type):
                diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts = 'yes'

            # 获取事件描述字符串
            for item in all_event_info_type.keys()[4:len(all_event_info_type.keys())]:
                line_18_items[line_18_items_dic[item]] = str(all_event_info_type.get_value(event_id, item))

            event_desc_str = '\t'.join(line_18_items)
            event_str = '\t'.join(
                [str(e) for e in [event_id, type, novel, old, all_event_info_type.get_value(event_id, 'GeneID'),
                                  all_event_info_type.get_value(event_id, 'geneSymbol'),
                                  all_event_info_type.get_value(event_id, 'chr'),
                                  all_event_info_type.get_value(event_id, 'strand'), event_desc_str]])
            if diff_ReadsOnTargetAndJunctionCounts == 'no':
                ReadsOnTargetAndJunctionCounts_str = '\t'.join(np.repeat('null', 17).tolist())
            else:
                ReadsOnTargetAndJunctionCounts_str = '\t'.join(
                    [str(e) for e in ReadsOnTargetAndJunctionCounts_info_type.ix[event_id]._values[
                                     -17:len(ReadsOnTargetAndJunctionCounts_info_type.keys())]])
            if diff_JunctionCountOnly == 'no':
                JunctionCountOnly_str = '\t'.join(np.repeat('null', 17).tolist())
            else:
                JunctionCountOnly_str = '\t'.join(
                    [str(e) for e in JunctionCountOnly_info_type.ix[event_id]._values[
                                     -17:len(JunctionCountOnly_info_type.keys())]])
            fw.write('\t'.join(
                [event_str, diff_JunctionCountOnly, JunctionCountOnly_str, diff_ReadsOnTargetAndJunctionCounts,
                 ReadsOnTargetAndJunctionCounts_str, diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts,
                 diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts]) + '\n')

    fw.close()


def format_line(infile, outfile):
    fw = open(outfile, 'w')
    for line in open(infile):
        pass

    fw.close()


def process_single_rmats_output_dir(root):
    big_detail_file = os.path.join(root, 'all_events_detail_big_table.txt')
    stat_psi = os.path.join(root, 'psi_stats.file.txt')
    events_stats_file = os.path.join(root, 'event_stats.file.txt')

    # 检查输出结果文件夹合理性
    alter_ref_dic = check_rmats_out_dir(root)  # alter_ref_dic

    # 修饰结果文件 在ID 上加上事件类型
    s1 = time.time()
    for in_file, out_file in alter_ref_dic.items():
        # format_line(infile=in_file, outfile=out_file[0])
        # modify_id_for_txt(infile=out_file[0], outfile=out_file[1])
        # modify_id_for_txt(infile=in_file, outfile=out_file[1])
        modify_id_for_txt(infile=in_file, outfile=out_file)
    s2 = time.time()
    print('转换ID任务完成，用时：{}'.format(s2 - s1))

    print('开始获取事件基本信息统计')
    event_info_dic = get_event_stats(files=alter_ref_dic.values())
    s3 = time.time()
    print('获取事件基本信息统计结束，用时：{}'.format(s3 - s2))
    print('开始写事件统计文件')
    write_event_stat_file(event_stat_f=events_stats_file, event_stat_dic=event_info_dic)
    s4 = time.time()
    print('写事件统计文件结束，用时：{}'.format(s4 - s3))
    print('开始读psi info dic')
    psi_info_dic = get_event_psi_dic(event_info_dic)
    s5 = time.time()
    print('读psi info dic结束，用时：{}'.format(s5 - s4))
    write_psi_detail_file(file_src_dic=event_info_dic)
    s6 = time.time()
    print('写psi细节文件结束，用时：{}'.format(s6 - s5))
    write_psi_stats_file(stat_file=stat_psi, psi_dic=psi_info_dic)
    s7 = time.time()
    print('写psi统计文件结束，用时：{}'.format(s7 - s6))
    write_big_detail_file(big_file=big_detail_file, out=root)
    s8 = time.time()
    print('写大文件结束，用时：{}'.format(s8 - s7))
    print('程序总用时：{}'.format(s8 - s1))


if __name__ == '__main__':
    process_single_rmats_output_dir(
        # root='/mnt/ilustre/users/sanger-dev/workspace/'
        #      '20170502/Single_rmats_module_linfang_new/Rmats/RmatsBam/output')
        root='/mnt/ilustre/users/sanger-dev/workspace/20170706/Rmats_demo1_1087_8595/RmatsBam/output')
