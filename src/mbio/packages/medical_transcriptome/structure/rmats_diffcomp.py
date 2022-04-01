# -*- coding: utf-8 -*-
# __author__ = 'gudeqing,qinjincheng'

from mbio.packages.medical_transcriptome.functions import pkgsfuncdeco
import os
from biocluster.config import Config
from bson.objectid import ObjectId
import pickle
import itertools
from collections import OrderedDict
import pandas as pd
import random

@pkgsfuncdeco
def main(args):
    info_dict = get_info_dict(args.input,args.version)
    file_tmp = os.path.join(os.getcwd(), 'diffcomp.tmp')
    if get_finaltable(args.input, info_dict, args.method, args.psi, args.cutoff, file_tmp):
        if get_output(file_tmp, info_dict, args.output):
            print 'succeed in exporting {}'.format(args.output)

@pkgsfuncdeco
def get_info_dict(paths_file,version):
    dct = {os.path.basename(line.strip()).rstrip('.txt'): dict() for line in open(paths_file)}
    db = Config().get_mongo_client(mtype='medical_transcriptome',db_version =version)[Config().get_mongo_dbname('medical_transcriptome',db_version =version)]
    for main_id in dct:
        dct[main_id]['compare_plan'] = db['sg_splicing_rmats'].find_one({'main_id': ObjectId(main_id)})['compare_plan']
        documents = list(db['sg_splicing_rmats_detail'].find({'splicing_id': ObjectId(main_id)}))
        for d in documents:
            dct[main_id].update(get_event_transcripts(d))
    else:
        return dct

def get_event_transcripts(document):
    k1, k2 = {'SE': ('InclusionTranscripts', 'SkippingTranscripts'),
              'RI': ('RetainTranscripts', 'AbandonTranscripts'),
              'MXE': ('1stExonTranscripts', '2ndExonTranscripts'),
              'A5SS': ('LongExonTranscripts', 'ShortExonTranscripts'),
              'A3SS': ('LongExonTranscripts', 'ShortExonTranscripts')}[document['type']]
    return {document['event_id']: (document[k1], document[k2])}

@pkgsfuncdeco
def get_finaltable(tablelist, infodict, significant_method, delta_psi, significant_cutoff, finaltable):
    title_detail = list()
    results = list()
    comparisons = list()
    for line in open(tablelist):
        result = line.strip()
        comparison = infodict[os.path.basename(result).rstrip('.txt')]['compare_plan']
        results.append(result)
        comparisons.append(comparison)
    uniq_event_new = diffgroup_stat(results, significant_method, delta_psi, significant_cutoff)
    detail_list = [
        '_AS_ID', '_diff_significant_JC', '_delta_PSI_JC', '_Pvalue_JC', '_FDR_JC',
        '_diff_significant_JCEC', '_delta_PSI_JCEC','_Pvalue_JCEC','_FDR_JCEC'
    ]
    header_list = list(itertools.product(comparisons, detail_list))
    for i in header_list:
        single = ''.join(list(i))
        title_detail.append(single)
    pickle.dump(uniq_event_new, open('uniq_event_new.pk', 'w'))
    with open(finaltable, 'w') as f:
        f.write('New_ID\tType\tGene_ID\tGene_name\t{}\n'.format('\t'.join(title_detail)))
        for k in uniq_event_new.keys():
            comp_result = list()
            for num, items in enumerate(results):
                if num in uniq_event_new[k][1]:
                    comp_result.insert(num, uniq_event_new[k].pop(2))
                else:
                    comp_result.insert(num, list('_________'))
            result_all = uniq_event_new[k][0:1] + list(k)[0:3] + list(itertools.chain(*comp_result))
            f.write('{}\n'.format('\t'.join(result_all)))
    return True

@pkgsfuncdeco
def diffgroup_stat(results, significant_method, delta_psi, significant_cutoff):
    uniq_event_ret = get_uniqid(results)
    for detail_file in results:
        print 'start processing {}'.format(detail_file)
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
                    if tuple(uniq_site) in uniq_event_ret.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(
                                significant_value) < significant_cutoff:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) >delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if line.split()[1] == 'A5SS':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                    if tuple(uniq_site) in uniq_event_ret.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(
                                significant_value) < significant_cutoff:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff :
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if line.split()[1] == 'MXE':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[22:26]
                    if tuple(uniq_site) in uniq_event_ret.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(
                                significant_value) < significant_cutoff:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if line.split()[1] == 'RI':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[20:22]
                    if tuple(uniq_site) in uniq_event_ret.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff :
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(
                                significant_value) < significant_cutoff:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff :
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if line.split()[1] == 'SE':
                    uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[8:14]
                    if tuple(uniq_site) in uniq_event_ret.keys():
                        if line.split()[39] == 'null':
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['null'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        elif abs(float(line.split()[39])) > delta_psi and float(significant_value) < significant_cutoff :
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['yes'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                        else:
                            if line.split()[57] == 'null':
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'null'] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[57])) > delta_psi and float(
                                    significant_value_all) < significant_cutoff:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'yes'] + line.split()[57:58] + line.split()[51:53])
                            else:
                                uniq_event_ret[tuple(uniq_site)].append(
                                    line.split()[0:1] + ['no'] + line.split()[39:40] + line.split()[33:35] + [
                                        'no'] + line.split()[57:58] + line.split()[51:53])
                if count % 1000 == 0:
                    print '{} lines have been inspected'.format(count)
                count += 1
            else:
                print 'succeed in processing {}'.format(detail_file)
    else:
        return uniq_event_ret

@pkgsfuncdeco
def get_uniqid(result_list):
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

@pkgsfuncdeco
def get_output(ifile, dct, ofile):
    rows = list()
    for n, row in pd.read_table(ifile).iterrows():
        if len({row['{}_AS_ID'.format(dct[main_id]['compare_plan'])] for main_id in dct}) == 1:
            splicing_id = dct.keys()[random.randint(0, len(dct) - 1)]
            event_id = row['{}_AS_ID'.format(dct[splicing_id]['compare_plan'])]
            row['Transcripts_1'], row['Transcripts_2'] = dct[splicing_id][event_id]
            rows.append(row.fillna('_'))
    else:
        pd.DataFrame(rows).to_csv(ofile, sep='\t', index=False)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Combine rMATS big tables')
    parser.add_argument('-i', '--input', dest='input', required=True, help='input rMATS big table paths list')
    parser.add_argument('-p', '--psi', dest='psi', required=True, type=float, help='idelta PSI')
    parser.add_argument('-m', '--method', dest='method', required=True, help='significant method')
    parser.add_argument('-c', '--cutoff', dest='cutoff', required=True, type=float, help='significant cutoff')
    parser.add_argument('-o', '--output', dest='output', required=True, help='output DIFFCOMP table')
    parser.add_argument('-v', '--version', dest='version', required=True, type=int, help='db_version')
    args = parser.parse_args()

    main(args)
