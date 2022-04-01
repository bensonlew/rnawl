# -*- coding: utf-8 -*-
# __author__ = 'jinlinfang,shicaiping,gudeqing,qinjincheng'

import re
import os
import Bio
import subprocess
from collections import defaultdict
import pandas as pd
import numpy as np
import glob
from mbio.packages.medical_transcriptome.functions import pkgsfuncdeco

global LEGAL_EVENT_TYPE
global PRIMARY_DIRS
global EVENT_FILES
global MATS_FILES

LEGAL_EVENT_TYPE = {'A3SS', 'A5SS', 'MXE', 'RI', 'SE'}
OUTPUT_FILES = {'fromGTF.A3SS.txt', 'fromGTF.A5SS.txt', 'fromGTF.MXE.txt',
                'fromGTF.novelEvents.A3SS.txt', 'fromGTF.novelEvents.A5SS.txt',
                'fromGTF.novelEvents.MXE.txt', 'fromGTF.novelEvents.RI.txt',
                'fromGTF.novelEvents.SE.txt', 'fromGTF.RI.txt', 'fromGTF.SE.txt',
                'A3SS.MATS.JC.txt', 'A3SS.MATS.JCEC.txt','A5SS.MATS.JC.txt',
                'A5SS.MATS.JCEC.txt', 'MXE.MATS.JC.txt','MXE.MATS.JCEC.txt',
                'RI.MATS.JC.txt', 'RI.MATS.JCEC.txt', 'SE.MATS.JC.txt', 'SE.MATS.JCEC.txt'}
EVENT_DESC_ITEMS = {'A3SS': ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
                    'A5SS': ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
                    'SE': ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
                    'MXE': ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES",
                            "upstreamEE", "downstreamES", "downstreamEE"],
                    'RI': ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES",
                           "downstreamEE"]}

@pkgsfuncdeco
def modify_id_for_txt(input_gtf, infile, outfile):
    '''
    this function is deprecated
    '''
    chr_set = [e.strip() for e in subprocess.check_output('awk -F \'\\t\'  \'$0!~/^#/{print $1}\' %s  | uniq | sort |uniq'
                                                          % input_gtf, shell=True).strip().split('\n')]
    print chr_set
    f_name_match = re.match(r'^.*?(A3SS|A5SS|MXE|RI|SE).*?\.txt$', os.path.basename(infile).strip())
    if f_name_match:
        tmp = outfile + '.temp'
        content = re.sub(r'\"', '', open(infile).read())
        open(tmp, 'w').write(content)
        event_type = f_name_match.group(1)
        cmd = """awk -F \'\\t\'  \'{if ($1~/[0-9]+/) printf  "%s"$0"\\n"; else print $0}\' %s > %s""" % (
            event_type + "_", tmp, outfile)
        print cmd
        subprocess.call(cmd, shell=True)
        with open (tmp, "r") as f: # 如果染色体编号不带chr，将改为与gtf文件染色体编号一致
            head = f.readline()
            a = f.readline()
            if a:
                chr = a.split("\t")[3]
                if chr not in chr_set:
                    cmd1 = """\n sed -i \'2,$s/\\tchr/\\t/\' %s""" %(outfile)
                    print cmd1
                    subprocess.call(cmd1, shell=True)
    else:
        raise Exception('输入的rmats结果文件名字有误：%s' % (infile))

@pkgsfuncdeco
def check_rmats_out_dir(dirpath):
    if OUTPUT_FILES.difference(set(os.listdir(dirpath))):
        raise Exception("不完整的rmats结果目录结构")
    else:
        # get event id list for filtering purpose. added by gdq
        event_types = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']
        event_id_dict = dict()
        for each in event_types:
            target_file = each + '.MATS.JCEC.txt'
            target_file = os.path.join(dirpath, target_file)
            with open(target_file) as f:
                f.readline()
                id_list = [int(x.split('\t')[0]) for x in f if x.strip()]
            event_id_dict[each] = id_list

        def drop_quotations(file):
            lines = [line.replace('"', '') for line in open(file)]
            open(file, 'w').writelines(lines)
        # modify the files of fromGTF*.txt
        for each in event_types:
            # print 'start check AS type -> {}'.format(each)
            target_file1 = 'fromGTF.{}.txt'.format(each)
            target_file2 = 'fromGTF.novelEvents.{}.txt'.format(each)
            target_file1 = os.path.join(dirpath, target_file1)
            target_file2 = os.path.join(dirpath, target_file2)
            id_list = event_id_dict[each]
            # modify 'fromGTF.{}.txt'
            drop_quotations(target_file1)
            target_file_table = pd.read_table(target_file1, index_col=0, header=0)
            target_file_table = target_file_table.query("index in @id_list")
            target_file_table.to_csv(target_file1, sep='\t')
            # modify fromGTF.novelEvents.{}.txt
            drop_quotations(target_file2)
            target_file2_table = pd.read_table(target_file2, index_col=0, header=0)
            target_file2_table = target_file2_table.query("index in @id_list")
            target_file2_table.to_csv(target_file2, sep='\t')

    to_be_altered_id_files = set([os.path.join(dirpath, name) for name in os.listdir(dirpath) if
                                  re.match(r'^fromGTF\.(novelEvents\.)?(A3SS|A5SS|MXE|RI|SE)\.txt$', name.strip())] +
                                 [os.path.join(dirpath, name) for name in os.listdir(dirpath) if
                                  re.match(r'^(A3SS|A5SS|MXE|RI|SE)\.MATS\.(JCEC|JC)\.txt$',
                                     name.strip())])
    info_dic = defaultdict(str)
    for f in to_be_altered_id_files:
        alter_id_file = os.path.join(os.path.dirname(f),
                                     re.match(r'^(\S+)\.txt$', os.path.basename(f).strip()).group(1) + '.alter_id.txt')
        info_dic[f] = alter_id_file
    else:
        return info_dic

@pkgsfuncdeco
def generate_event_id(infile, outfile):
    required_basenames = ['A3SS.MATS.JC.txt', 'A3SS.MATS.JCEC.txt',
                          'A5SS.MATS.JC.txt', 'A5SS.MATS.JCEC.txt',
                          'MXE.MATS.JC.txt', 'MXE.MATS.JCEC.txt',
                          'RI.MATS.JC.txt', 'RI.MATS.JCEC.txt',
                          'SE.MATS.JC.txt', 'SE.MATS.JCEC.txt',
                          'fromGTF.A3SS.txt', 'fromGTF.A5SS.txt',
                          'fromGTF.MXE.txt', 'fromGTF.RI.txt',
                          'fromGTF.SE.txt', 'fromGTF.novelEvents.A3SS.txt',
                          'fromGTF.novelEvents.A5SS.txt', 'fromGTF.novelEvents.MXE.txt',
                          'fromGTF.novelEvents.RI.txt', 'fromGTF.novelEvents.SE.txt']
    assert os.path.basename(infile) in required_basenames
    print 'INFO: start reading {}'.format(infile)
    event_type = re.match(r'^.*?(A3SS|A5SS|MXE|RI|SE).*?\.txt$', os.path.basename(infile)).group(1)
    def drop_quotations(file):
        lines = [line.replace('"', '') for line in open(file)]
        open(file, 'w').writelines(lines)
    drop_quotations(infile)
    df = pd.read_table(infile)
    df['ID'] = df['ID'].apply(lambda x: '{}_{}'.format(event_type, x))
    df['chr'] = df['chr'].apply(lambda x: x.lstrip('chr'))
    df = df.fillna('NA')
    df.to_csv(outfile, sep='\t', index=False)

def isnumber(aString):
    try:
        float(aString)
        return True
    except:
        return False


def union(set_list):
    ret = set()
    ret.update(*set_list)
    return ret

@pkgsfuncdeco
def get_event_stats(files, pvalue_fdr='fdr', fdr=0.05, psi=0):
    dct = dict.fromkeys(LEGAL_EVENT_TYPE)
    for i in LEGAL_EVENT_TYPE:
        dct[i] = {'JunctionCountOnly_event_id_set': set(), 'ReadsOnTargetAndJunctionCounts_event_id_set': set(),
                  'JunctionCountOnly_event_id_set_no': 0, 'ReadsOnTargetAndJunctionCounts_event_id_set_no': 0,
                  'JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set': set(),
                  'JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no': 0,
                  'JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set': set(),
                  'JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no': 0}
    for file in files:
        m = re.match(r'(\S+?)\.MATS.(JC|JCEC)\.alter_id\.txt', os.path.basename(file))
        if m:
            print 'start reading {}'.format(file)
            event_type = m.group(1)
            dct[event_type]['{}_file'.format(m.group(2))] = file
            df = pd.read_table(file)
            if pvalue_fdr.lower() == 'fdr':
                fldf = df[df['FDR'] <= fdr]
                fldf = fldf[abs(fldf['IncLevelDifference']) >= psi]
            else:
                fldf = df[df['PValue']<=fdr]
                fldf = fldf[abs(fldf['IncLevelDifference']) >= psi]
            if m.group(2) == 'JC':
                dct[event_type]['JunctionCountOnly_event_id_set_no'] = len(set(fldf['ID']))
                dct[event_type]['JunctionCountOnly_event_id_set'] = set(fldf['ID'])
            elif m.group(2) == 'JCEC':
                dct[event_type]['ReadsOnTargetAndJunctionCounts_event_id_set_no'] = len(set(fldf['ID']))
                dct[event_type]['ReadsOnTargetAndJunctionCounts_event_id_set'] = set(fldf['ID'])
    for j in LEGAL_EVENT_TYPE:
        dct[j]['JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set'] = \
            dct[j]['JunctionCountOnly_event_id_set'] & dct[j]['ReadsOnTargetAndJunctionCounts_event_id_set']
        dct[j]['JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no'] = \
            len(dct[j]['JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set'])
        dct[j]['JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set'] = \
            dct[j]['JunctionCountOnly_event_id_set'] | dct[j]['ReadsOnTargetAndJunctionCounts_event_id_set']
        dct[j]['JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no'] = \
            len(dct[j]['JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set'])
    dct['total_JunctionCountOnly_event_id_set_no'] = sum(
        [dct[k]['JunctionCountOnly_event_id_set_no'] for k in LEGAL_EVENT_TYPE])
    dct['total_JunctionCountOnly_event_id_set'] = union(
        [dct[k]['JunctionCountOnly_event_id_set'] for k in LEGAL_EVENT_TYPE])
    dct['total_ReadsOnTargetAndJunctionCounts_event_id_set_no'] = sum(
        [dct[k]['ReadsOnTargetAndJunctionCounts_event_id_set_no'] for k in LEGAL_EVENT_TYPE])
    dct['total_ReadsOnTargetAndJunctionCounts_event_id_set'] = union(
        [dct[k]['ReadsOnTargetAndJunctionCounts_event_id_set'] for k in LEGAL_EVENT_TYPE])
    return dct

@pkgsfuncdeco
def get_event_type(files):
    dct = dict.fromkeys(LEGAL_EVENT_TYPE)
    for i in LEGAL_EVENT_TYPE:
        dct[i] = {'all_event_id_no': 0, 'all_event_id_set': set(),
                  'old_event_id_no': 0, 'old_event_id_set': set(),
                  'novel_event_id_no': 0, 'novel_event_id_set': set()}
    for file in files:
        m = re.match(r'fromGTF\.(novelEvents\.)?(A3SS|A5SS|MXE|RI|SE)\.alter_id\.txt', os.path.basename(file))
        if m:
            print 'INFO: start reading {}'.format(file)
            event_type = m.group(2)
            df = pd.read_table(file)
            if m.group(1):
                dct[event_type]['novel_event_file'] = file
                dct[event_type]['novel_event_id_no'] = len(set(df['ID']))
                dct[event_type]['novel_event_id_set'] = set(df['ID'])
            else:
                dct[event_type]['all_event_file'] = file
                dct[event_type]['all_event_id_no'] = len(set(df['ID']))
                dct[event_type]['all_event_id_set'] = set(df['ID'])
    for j in LEGAL_EVENT_TYPE:
        dct[j]['old_event_id_set'] = dct[j]['all_event_id_set'] - dct[j]['novel_event_id_set']
        dct[j]['old_event_id_set_no'] = len(dct[j]['old_event_id_set'])
    dct['total_as_events_no'] = sum([dct[k]['all_event_id_no'] for k in LEGAL_EVENT_TYPE])
    dct['total_as_events_id_set'] = union([dct[k]['all_event_id_set'] for k in LEGAL_EVENT_TYPE])
    dct['total_as_novel_events_no'] = sum([dct[k]['novel_event_id_no'] for k in LEGAL_EVENT_TYPE])
    dct['total_as_novel_events_id_set'] = union([dct[k]['novel_event_id_set'] for k in LEGAL_EVENT_TYPE])
    return dct

@pkgsfuncdeco
def get_event_psi_dic(dct, pvalue_fdr='fdr', fdr=0.05, psi=0):
    psi_dic = dict.fromkeys(LEGAL_EVENT_TYPE)
    psi_dic['SAMPLE_1'] = {'JunctionCountOnly': {'exclusion_total': 0, 'inclusion_total': 0, 'total': 0},
                           'ReadsOnTargetAndJunctionCounts': {'exclusion_total': 0, 'inclusion_total': 0, 'total': 0}}
    psi_dic['SAMPLE_2'] = {'JunctionCountOnly': {'exclusion_total': 0, 'inclusion_total': 0, 'total': 0},
                           'ReadsOnTargetAndJunctionCounts': {'exclusion_total': 0, 'inclusion_total': 0, 'total': 0}}
    psi_dic['total'] = {'JunctionCountOnly': 0, 'ReadsOnTargetAndJunctionCounts': 0}
    for k in LEGAL_EVENT_TYPE:
        psi_dic[k] = {'JunctionCountOnly_data': {},
                      'ReadsOnTargetAndJunctionCounts_data': {},
                      'SAMPLE_1': {'JunctionCountOnly': {'exclusion': 0, 'inclusion': 0, 'total': 0},
                                   'ReadsOnTargetAndJunctionCounts': {'exclusion': 0, 'inclusion': 0, 'total': 0}},
                      'SAMPLE_2': {'JunctionCountOnly': {'exclusion': 0, 'inclusion': 0, 'total': 0},
                                   'ReadsOnTargetAndJunctionCounts': {'exclusion': 0, 'inclusion': 0, 'total': 0}}}
    for event_type in LEGAL_EVENT_TYPE:
        print 'start processing event type -> ({})'.format(event_type)
        print 'start reading {}'.format(dct[event_type]['JC_file'])
        print 'start reading {}'.format(dct[event_type]['JCEC_file'])
        JC_file = dct[event_type]['JC_file']
        JCEC_file = dct[event_type]['JCEC_file']
        psi_dic['JC_file'] = JC_file
        psi_dic['JCEC_file'] = JCEC_file

        jc_data = pd.read_table(JC_file, index_col=0, sep='\t', dtype={
            'Pvalue': np.float64, 'FDR': np.float64, 'IncLevelDifference': np.float64
        })
        tmp_jc = jc_data[['PValue', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference']]
        if pvalue_fdr.lower() == 'fdr':
            df_filter = tmp_jc[tmp_jc['FDR']<=fdr]
            jc_df = df_filter[abs(df_filter['IncLevelDifference']) >= psi]
        else:
            df_filter = tmp_jc[tmp_jc['PValue']<=fdr]
            jc_df = df_filter[abs(df_filter['IncLevelDifference']) >= psi]
        jc_it = jc_df.iterrows()
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
        all_data = pd.read_table(JCEC_file, index_col=0, sep='\t',  dtype={
            'PValue': np.float64, 'FDR': np.float64, 'IncLevelDifference': np.float64
        })
        tmp_all = all_data[['PValue', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference']]
        all_df = tmp_all[tmp_all['FDR']<=0.05]
        all_it = all_df.iterrows()
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

@pkgsfuncdeco
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
    data = pd.read_table(mats_file, sep='\t', dtype={
        'PValue': float, 'FDR': float, 'IncLevelDifference': float, 'IncFormLen': int, 'SkipFormLen': int, 'ID.1': str
    })
    fw = open(new_file, 'wb')
    fw.write('{}\t{}\n'.format('\t'.join(data.keys()[0:len(data.keys()) - 1]), '\t'.join(
        ['average_IncLevel1', 'average_IncLevel2', 'IncLevelDifference', 'increase_inclusion_SAMPLE1',
         'increase_exclusion_SAMPLE1', 'increase_inclusion_SAMPLE2', 'increase_exclusion_SAMPLE2'])))
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

@pkgsfuncdeco
def write_psi_detail_file(file_src_dic):
    for as_type in LEGAL_EVENT_TYPE:
        JunctionCountOnly_file = file_src_dic[as_type]['JC_file']
        ReadsOnTargetAndJunctionCounts_file = file_src_dic[as_type]['JCEC_file']
        new_JunctionCountOnly_file = os.path.join(
            os.path.dirname(JunctionCountOnly_file),
            re.sub(r'(\S+)\.txt', '\g<1>.psi_info.txt', os.path.basename(JunctionCountOnly_file))
        )
        new_ReadsOnTargetAndJunctionCounts_file = os.path.join(
            os.path.dirname(ReadsOnTargetAndJunctionCounts_file),
            re.sub(r'(\S+)\.txt', '\g<1>.psi_info.txt', os.path.basename(ReadsOnTargetAndJunctionCounts_file))
        )
        add_psi_info(JunctionCountOnly_file, new_file=new_JunctionCountOnly_file)
        add_psi_info(ReadsOnTargetAndJunctionCounts_file, new_file=new_ReadsOnTargetAndJunctionCounts_file)
        file_src_dic[as_type]['JunctionCountOnly_add_psi_file'] = new_JunctionCountOnly_file
        file_src_dic[as_type]['ReadsOnTargetAndJunctionCounts_add_psi_file'] = new_ReadsOnTargetAndJunctionCounts_file
    return file_src_dic

@pkgsfuncdeco
def write_event_stat_file(event_stat_dic, event_stat_f):
    fw = open(event_stat_f, 'wb')
    fw.write('stat_item\tvalue\n')
    for info_key in event_stat_dic.keys():
        if info_key.strip() in LEGAL_EVENT_TYPE:
            for type_info_key in event_stat_dic[info_key].keys():
                if re.match(r'.+_no$', type_info_key.strip()):
                    fw.write('{}\t{}\n'.format(
                        '{}_{}'.format(info_key, type_info_key), event_stat_dic[info_key][type_info_key]
                    ))
                    continue
        if re.match(r'.+_no$', info_key.strip()):
            fw.write('{}\t{}\n'.format(info_key, event_stat_dic[info_key]))
            continue
    fw.close()

@pkgsfuncdeco
def write_event_type_file(event_stat_dic, event_stat_f):
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

@pkgsfuncdeco
def write_big_detail_file(out, big_file, pvalue_fdr="fdr", fdr=0.05, psi=0):
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
        all_event_file_type = os.path.join(out, 'fromGTF.%s.alter_id.txt' % (type))
        novel_event_file_type = os.path.join(out, 'fromGTF.novelEvents.%s.alter_id.txt' % (type))
        JunctionCountOnly_file_type = os.path.join(out,
                                                   '%s.MATS.JC.alter_id.psi_info.txt' % (
                                                       type))
        ReadsOnTargetAndJunctionCounts_file_type = os.path.join(out,
                                                                '%s.MATS.JCEC.alter_id.psi_info.txt' % (
                                                                    type))
        all_event_info_type = pd.read_table(all_event_file_type, index_col=0, sep='\t')
        all_event_set_type = set(all_event_info_type.index.tolist())
        novel_event_set_type = set(
            pd.read_table(novel_event_file_type, index_col=0, sep='\t').index.tolist())
        JunctionCountOnly_info_type = pd.read_table(JunctionCountOnly_file_type, sep='\t',
                                                        dtype={'FDR': str, 'PValue': str, 'IncLevelDifference': str,
                                                               'average_IncLevel2': str, 'average_IncLevel1': str},
                                                        index_col=0)
        JunctionCountOnly_event_set_type = set(JunctionCountOnly_info_type.index.tolist())
        ReadsOnTargetAndJunctionCounts_info_type = pd.read_table(ReadsOnTargetAndJunctionCounts_file_type,
                                                                     dtype={'FDR': str, 'PValue': str,
                                                                            'IncLevelDifference': str,
                                                                            'average_IncLevel2': str,
                                                                            'average_IncLevel1': str}, sep='\t',
                                                                     index_col=0)
        ReadsOnTargetAndJunctionCounts_event_set_type = set(ReadsOnTargetAndJunctionCounts_info_type.index.tolist())

        for event_id in all_event_set_type:
            novel = 'no'
            old = 'yes'
            diff_ReadsOnTargetAndJunctionCounts = 'null'
            diff_JunctionCountOnly = 'null'
            diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts = 'null'
            diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts = 'null'

            if event_id in novel_event_set_type:
                novel = 'yes'
                old = 'no'
            if event_id in ReadsOnTargetAndJunctionCounts_event_set_type:
                if pvalue_fdr.lower() == 'fdr':
                    if float(ReadsOnTargetAndJunctionCounts_info_type.loc[event_id]['FDR']) <= fdr and abs(float(ReadsOnTargetAndJunctionCounts_info_type.loc[event_id]['IncLevelDifference'])) >= psi:
                        diff_ReadsOnTargetAndJunctionCounts = 'yes'
                    else:
                        diff_ReadsOnTargetAndJunctionCounts = 'no'
                else:
                    if float(ReadsOnTargetAndJunctionCounts_info_type.loc[event_id]['PValue']) <= fdr and abs(float(ReadsOnTargetAndJunctionCounts_info_type.loc[event_id]['IncLevelDifference'])) >= psi:
                        diff_ReadsOnTargetAndJunctionCounts = 'yes'
                    else:
                        diff_ReadsOnTargetAndJunctionCounts = 'no'
            if event_id in JunctionCountOnly_event_set_type:
                if pvalue_fdr.lower() == 'fdr':
                    if float(JunctionCountOnly_info_type.loc[event_id]['FDR']) <= fdr and abs(float(JunctionCountOnly_info_type.loc[event_id]['IncLevelDifference'])) >= psi:
                        diff_JunctionCountOnly = 'yes'
                    else:
                        diff_JunctionCountOnly = 'no'
                else:
                    if float(JunctionCountOnly_info_type.loc[event_id]['PValue']) <= fdr and abs(float(JunctionCountOnly_info_type.loc[event_id]['IncLevelDifference'])) >= psi:
                        diff_JunctionCountOnly = 'yes'
                    else:
                        diff_JunctionCountOnly = 'no'
            if diff_JunctionCountOnly == 'yes' and diff_ReadsOnTargetAndJunctionCounts == 'yes':
                diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts = 'yes'
                diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts = 'yes'
            elif diff_JunctionCountOnly == 'yes' or diff_ReadsOnTargetAndJunctionCounts == 'yes':
                diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts = 'yes'
                diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts = 'no'
            else:
                diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts = 'no'
                diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts = 'no'
            '''
            if (event_id in JunctionCountOnly_event_set_type) and (
                        event_id in ReadsOnTargetAndJunctionCounts_event_set_type):
                diff_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts = 'yes'
            if (event_id in JunctionCountOnly_event_set_type) or (
                        event_id in ReadsOnTargetAndJunctionCounts_event_set_type):
                diff_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts = 'yes'
            '''

            # 获取事件描述字符串
            for item in all_event_info_type.keys()[4:len(all_event_info_type.keys())]:
                line_18_items[line_18_items_dic[item]] = str(all_event_info_type.at[event_id, item])

            event_desc_str = '\t'.join(line_18_items)
            event_str = '\t'.join(
                [str(e) for e in [event_id, type, novel, old, all_event_info_type.at[event_id, 'GeneID'],
                                  all_event_info_type.at[event_id, 'geneSymbol'],
                                  all_event_info_type.at[event_id, 'chr'],
                                  all_event_info_type.at[event_id, 'strand'], event_desc_str]])
            if diff_ReadsOnTargetAndJunctionCounts == 'null':
                ReadsOnTargetAndJunctionCounts_str = '\t'.join(np.repeat('null', 17).tolist())
            else:
                ReadsOnTargetAndJunctionCounts_str = '\t'.join(
                    [str(e) for e in ReadsOnTargetAndJunctionCounts_info_type.ix[event_id]._values[
                                     -17:len(ReadsOnTargetAndJunctionCounts_info_type.keys())]])
            if diff_JunctionCountOnly == 'null':
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

@pkgsfuncdeco
def process_single_rmats_output_dir(root, input_gtf, pvalue_fdr='fdr', fdr=0.05, psi=0):
    # 规定文件位置以及定义变量
    all_table_file = os.path.join(root, 'all_events_detail_big_table.txt')
    psi_stats_file = os.path.join(root, 'psi_stats.file.txt')
    event_stats_file = os.path.join(root, 'event_stats.file.txt')
    event_type_file = os.path.join(root, 'event_type.file.txt')

    # 检查输出结果文件夹合理性
    alter_ref_dic = check_rmats_out_dir(root)

    # 修饰结果文件以及在ID前加上事件类型
    for infile, outfile in alter_ref_dic.items():
        generate_event_id(infile, outfile)

    # 获取事件基本信息统计
    event_info_dic = get_event_stats(files=alter_ref_dic.values(), pvalue_fdr=pvalue_fdr, fdr=fdr, psi=psi)
    event_type_dic = get_event_type(files=alter_ref_dic.values())

    # 导出事件统计文件
    write_event_stat_file(event_stat_f=event_stats_file, event_stat_dic=event_info_dic)
    write_event_type_file(event_stat_f=event_type_file, event_stat_dic=event_type_dic)

    # 获取PSI信息字典
    psi_info_dic = get_event_psi_dic(event_info_dic, pvalue_fdr=pvalue_fdr, fdr=fdr, psi=psi)

    # 导出PSI细节文件
    write_psi_detail_file(file_src_dic=event_info_dic)

    # 导出PSI统计文件
    write_psi_stats_file(stat_file=psi_stats_file, psi_dic=psi_info_dic)

    # 导出大文件
    write_big_detail_file(big_file=all_table_file, out=root, pvalue_fdr=pvalue_fdr, fdr=fdr, psi=psi)

def process_rmats_stat(root, pvalue_fdr, fdr, psi):
    files = glob.glob(os.path.join(root, '*.txt'))
    events_stats_file = os.path.join(root, 'event_stats.file.txt')
    psi_stats_file = os.path.join(root, 'psi_stats.file.txt')

    # 获取事件基本信息统计
    event_info_dic = get_event_stats(files=files, pvalue_fdr=pvalue_fdr, fdr=fdr, psi=psi)

    # 导出事件统计文件
    write_event_stat_file(event_stat_f=events_stats_file, event_stat_dic=event_info_dic)

    # 获取PSI信息字典
    psi_info_dic = get_event_psi_dic(event_info_dic, pvalue_fdr=pvalue_fdr, fdr=fdr, psi=psi)

    # 导出PSI细节文件
    write_psi_detail_file(file_src_dic=event_info_dic)

    # 导出PSI统计文件
    write_psi_stats_file(stat_file=psi_stats_file, psi_dic=psi_info_dic)

if __name__ == '__main__':
    import sys
    import argparse
    import unittest

    parser = argparse.ArgumentParser(description='Process rMATS last output')
    subparsers = parser.add_subparsers(help='sub-commands')

    parser_alter = subparsers.add_parser('alter', help='process rMATS raw output')
    parser_alter.add_argument('--root', action='store', required=True, dest='root',
                              help='raw output folder of rMATS')
    parser_alter.add_argument('--gtf', action='store', required=True, dest='gtf',
                              help='annotation file of genes and transcripts in GTF format')
    parser_alter.add_argument('--sig', action='store', choices=['fdr', 'pvalue'], required=True, dest='sig',
                              help='measure of significance testing')
    parser_alter.add_argument('--cutoff', action='store', type=float, required=True, dest='cutoff',
                              help='cutoff value of significance testing')
    parser_alter.add_argument('--psi', action='store', type=float, required=True, dest='psi',
                              help='cutoff value of splicing difference')

    parser_stats = subparsers.add_parser('stats', help='stats rMATS last output')
    parser_stats.add_argument('--root', action='store', required=True, dest='root',
                              help='raw output folder of rMATS')
    parser_stats.add_argument('--sig', action='store', choices=['fdr', 'pvalue'], required=True, dest='sig',
                              help='measure of significance testing')
    parser_stats.add_argument('--cutoff', action='store', type=float, required=True, dest='cutoff',
                              help='cutoff value of significance testing')
    parser_stats.add_argument('--psi', action='store', type=float, required=True, dest='psi',
                              help='cutoff value of splicing difference')

    args = parser.parse_args()

    class TestFunction(unittest.TestCase):
        '''
        This is test for the function. Just run script to do test.
        '''
        def test_alter(self):
            process_single_rmats_output_dir(
                root=args.root,
                input_gtf=args.gtf,
                pvalue_fdr=args.sig,
                fdr=args.cutoff,
                psi=args.psi
            )
        def test_stats(self):
            process_rmats_stat(
                root=args.root,
                input_gtf=args.gtf,
                pvalue_fdr=args.sig,
                fdr=args.cutoff,
                psi=args.psi
            )

    suite = unittest.TestSuite()
    if sys.argv[1] == 'alter':
        suite.addTests([TestFunction('test_alter')])
    elif sys.argv[1] == 'stats':
        suite.addTests([TestFunction('test_stats')])
    unittest.TextTestRunner(verbosity=2).run(suite)
