# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from optparse import OptionParser
import os
import sys
import pickle
from collections import OrderedDict
import itertools
import pandas as pd

parser = OptionParser(description='Combine rMATS output detail files')
parser.add_option('-i', '--input', dest='input', help='input LOC2NAME file')
parser.add_option('-m', '--method', dest='method', help='significant method')
parser.add_option('-p', '--psi', dest='psi', type=float, help='idelta PSI')
parser.add_option('-c', '--cutoff', dest='cutoff', type=float, help='significant cutoff')
parser.add_option('-t', '--type', dest='type', type=str, help='input TYPE file')
parser.add_option('-o', '--output', dest='output', help='output DIFFCOMP file')
(opts, args) = parser.parse_args()

def main(file_in, method, psi, cutoff, type_file, file_out):
    print 'INFO: start reading {}'.format(file_in)
    file_tmp = os.path.join(os.getcwd(), 'diffcomp.tmp')
    ret = get_finaltable(file_in, method, psi, cutoff, file_tmp)
    if ret:
        open(file_out, 'w').writelines([line.replace('null', '_') for line in open(file_tmp)])
        df = pd.read_table(file_out).merge(
            pd.read_table(type_file, names=['Gene_ID', 'rna_type'], usecols=[0, 2]), on='Gene_ID', how='left'
        )
        df.to_csv(file_out, sep='\t', index=None)
        if os.path.getsize(file_out) > 0:
            print 'INFO: succeed in exporting {}'.format(file_out)

def get_finaltable(loc2name, significant_method, delta_psi, significant_cutoff, finaltable):
    print 'INFO: start calling {}'.format(sys._getframe().f_code.co_name)
    title_detail = list()
    result_list = list()
    diffgroup = list()
    for line in open(loc2name):
        loc, name = line.strip().split('\t')
        result_list.append(loc)
        diffgroup.append(name)
    uniq_event_new = diffgroup_stat(result_list, significant_method, delta_psi, significant_cutoff)
    if uniq_event_new:
        print 'INFO: succeed in calling diffgroup_stat'
    detail_list = [
        '_AS_ID', '_diff_significant_JC', '_delta_PSI_JC', '_Pvalue_JC', '_FDR_JC',
        '_diff_significant_JCEC', '_delta_PSI_JCEC','_Pvalue_JCEC','_FDR_JCEC'
    ]
    header_list = list(itertools.product(diffgroup, detail_list))
    for i in header_list:
        single = ''.join(list(i))
        title_detail.append(single)
    pickle.dump(uniq_event_new, open('uniq_event_new.pk', 'w'))
    with open(finaltable, 'w') as f:
        f.write('New_ID\tType\tGene_ID\tGene_name\t{}\n'.format('\t'.join(title_detail)))
        for v in uniq_event_new.keys():
            comp_result = list()
            for num, items in enumerate(result_list):
                if num in uniq_event_new[v][1]:
                    comp_result.insert(num, uniq_event_new[v].pop(2))
                else:
                    comp_result.insert(num, list('_________'))
            result_all = uniq_event_new[v][0:1] + list(v)[0:3] + list(itertools.chain(*comp_result))
            f.write('{}\n'.format('\t'.join(result_all)))
    return True

def diffgroup_stat(result_list, significant_method, delta_psi, significant_cutoff):
    print 'INFO: start calling {}'.format(sys._getframe().f_code.co_name)
    uniq_event_new = get_uniqid(result_list)
    if uniq_event_new:
        print 'INFO: succeed in calling get_uniqid'
    for detail_file in result_list:
        print 'INFO: start processing {}'.format(detail_file)
        with open(detail_file) as f:
            _ = f.readline()
            count = 1
            for line in f:
                if significant_method.lower() == 'fdr':
                    significant_value = line.split()[34]
                    significant_value_all = line.split()[52]
                elif significant_method.lower() == 'pvalue':
                    significant_value = line.split()[33]
                    significant_value_all = line.split()[51]
                if line.split()[1] == 'A3SS':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                    if tuple(uniq_site) in uniq_event_new.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(
                                significant_value) < significant_cutoff:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) >delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if line.split()[1] == 'A5SS':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                    if tuple(uniq_site) in uniq_event_new.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(
                                significant_value) < significant_cutoff:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff :
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if line.split()[1] == 'MXE':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[22:26]
                    if tuple(uniq_site) in uniq_event_new.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(
                                significant_value) < significant_cutoff:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if line.split()[1] == 'RI':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[20:22]
                    if tuple(uniq_site) in uniq_event_new.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff :
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(
                                significant_value) < significant_cutoff:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff :
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if line.split()[1] == 'SE':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[8:14]
                    if tuple(uniq_site) in uniq_event_new.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(significant_value) < significant_cutoff :
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_new[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if count % 10000 == 0:
                    print 'INFO: {} lines have been inspected'.format(count)
                count += 1
    return uniq_event_new

def get_uniqid(result_list):
    print 'INFO: start calling {}'.format(sys._getframe().f_code.co_name)
    count = 0
    uniq_event = OrderedDict()
    for g_num, items in enumerate(result_list):
        with open(items) as f:
            _ = f.readline()
            for (num,line) in enumerate(f):
                if line.split()[1] == 'A3SS':
                    items_list = []
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                    if not uniq_event.has_key(tuple(uniq_site)):
                        uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                        items_list.append(g_num)
                        uniq_event[tuple(uniq_site)].append(items_list)
                    else:
                        uniq_event[tuple(uniq_site)][1].append(g_num)
                if line.split()[1] == 'A5SS':
                    items_list = []
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                    if not uniq_event.has_key(tuple(uniq_site)):
                        uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                        items_list.append(g_num)
                        uniq_event[tuple(uniq_site)].append(items_list)
                    else:
                        uniq_event[tuple(uniq_site)][1].append(g_num)
                if line.split()[1] == 'MXE':
                    items_list = []
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[22:26]
                    if not uniq_event.has_key(tuple(uniq_site)):
                        uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                        items_list.append(g_num)
                        uniq_event[tuple(uniq_site)].append(items_list)
                    else:
                        uniq_event[tuple(uniq_site)][1].append(g_num)
                if line.split()[1] == 'RI':
                    items_list = []
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[20:22]
                    if not uniq_event.has_key(tuple(uniq_site)):
                        uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                        items_list.append(g_num)
                        uniq_event[tuple(uniq_site)].append(items_list)
                    else:
                        uniq_event[tuple(uniq_site)][1].append(g_num)
                if line.split()[1] == 'SE':
                    items_list = []
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[8:14]
                    if not uniq_event.has_key(tuple(uniq_site)):
                        uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                        items_list.append(g_num)
                        uniq_event[tuple(uniq_site)].append(items_list)
                    else:
                        uniq_event[tuple(uniq_site)][1].append(g_num)
        count = count + num + 1
    return uniq_event

if __name__ == '__main__':
    if opts.input and opts.method and hasattr(opts, 'psi') and opts.cutoff and opts.type and opts.output:
        main(opts.input, opts.method, opts.psi, opts.cutoff, opts.type, opts.output)
    else:
        parser.print_help()