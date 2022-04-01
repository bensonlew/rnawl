# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import pickle

import numpy as np
import pandas as pd
from biocluster.config import Config
from bson.objectid import ObjectId

from mbio.packages.ref_rna_v3.functions import pkgsfuncdeco

# database = Config().get_mongo_client(mtype='medical_transcriptome')[Config().get_mongo_dbname('medical_transcriptome')]


def main(args):
    dct = {et: {'JC': {'s1in': set(), 's1ex': set()}, 'JCEC': {'s1in': set(), 's1ex': set()}} for et in
           ('SE', 'MXE', 'A3SS', 'A5SS', 'RI')}
    if args.version == 1:
        database = Config().get_mongo_client(mtype='medical_transcriptome',db_version =args.version)[Config().get_mongo_dbname('medical_transcriptome',db_version =args.version)]
    else:
        database = Config().get_mongo_client(mtype='medical_transcriptome')[Config().get_mongo_dbname('medical_transcriptome')]
    params = database['sg_splicing_rmats_stats'].find_one({'main_id': ObjectId(args.main_id)})['params']
    splicing_id = ObjectId(json.loads(params)['splicing_id'])
    cursor = database['sg_splicing_rmats_detail'].find({'splicing_id': splicing_id})

    detail_data = list()
    for document in cursor:
        detail_data.append([document['event_id'], document['inc_diff_jc'], document['inc_diff_all']])
    detail_df = pd.DataFrame(detail_data, columns=['AS ID', 'IncLevelDiffJC', 'IncLevelDiffJCEC'])
    detail_df = detail_df.set_index('AS ID')

    pickle.dump(detail_df, open('detail.df.pkl', 'wb'))
    detail_df['IncLevelDiffJC'] = detail_df['IncLevelDiffJC'].apply(lambda x: 0 if not isinstance(x, float) else x)
    detail_df['IncLevelDiffJCEC'] = detail_df['IncLevelDiffJCEC'].apply(lambda x: 0 if not isinstance(x, float) else x)
    detail_df['IncLevelDiffJC'] = detail_df['IncLevelDiffJC'].astype(float)
    detail_df['IncLevelDiffJCEC'] = detail_df['IncLevelDiffJCEC'].astype(float)

    for et in ('SE', 'MXE', 'A3SS', 'A5SS', 'RI'):
        for jt in ('JC', 'JCEC'):
            df = pd.read_table(os.path.join(args.root, '{}/{}.detail.xls'.format(jt, et)))
            sigcol = {
                'pvalue': {'JC': 'P Value JunctionCountOnly', 'JCEC': 'P Value ReadsOnTargetAndJunctionCounts'},
                'fdr': {'JC': 'FDR JunctionCountOnly', 'JCEC': 'FDR ReadsOnTargetAndJunctionCounts'}
            }[args.pvalue_fdr][jt]
            df = df[df[sigcol] <= args.fdr]

            df = df.set_index('AS ID')
            df = df.join(detail_df)
            df = df.reset_index()
            psicol = 'IncLevelDiff{}'.format(jt)

            pickle.dump(df, open('df.{}.{}.pkl'.format(jt, et), 'wb'))
            df = df[abs(df[psicol]) >= args.psi]
            for n, row in df.iterrows():
                if row[psicol] >= 0:
                    dct[et][jt]['s1in'].add(row['AS ID'])
                else:
                    dct[et][jt]['s1ex'].add(row['AS ID'])
    else:
        export_event_stats(dct, os.path.join(args.output, 'event_stats.file.txt'))
        export_psi_stats(dct, os.path.join(args.output, 'psi_stats.file.txt'))


@pkgsfuncdeco
def export_event_stats(dct, ofile):
    lines = ['stat_item\tvalue\n']
    total_jc_event_ids_no = 0
    total_jcec_event_ids_no = 0
    for et, jt2s1dct in dct.items():
        jc_event_ids = set()
        jc_event_ids.update(jt2s1dct['JC']['s1in'], jt2s1dct['JC']['s1ex'])
        jcec_event_ids = set()
        jcec_event_ids.update(jt2s1dct['JCEC']['s1in'], jt2s1dct['JCEC']['s1ex'])
        its_event_ids = np.intersect1d(list(jc_event_ids), list(jcec_event_ids))
        all_event_ids = set()
        all_event_ids.update(jc_event_ids, jcec_event_ids)
        lines.extend([
            '{}_JunctionCountOnly_event_id_set_no\t{}\n'.format(et, len(jc_event_ids)),
            '{}_ReadsOnTargetAndJunctionCounts_event_id_set_no\t{}\n'.format(et, len(jcec_event_ids)),
            '{}_JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no\t{}\n'.format(et, len(its_event_ids)),
            '{}_JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no\t{}\n'.format(et, len(all_event_ids))
        ])
        total_jc_event_ids_no += len(jc_event_ids)
        total_jcec_event_ids_no += len(jcec_event_ids)
    else:
        lines.extend([
            'total_JunctionCountOnly_event_id_set_no\t{}\n'.format(total_jc_event_ids_no),
            'total_ReadsOnTargetAndJunctionCounts_event_id_set_no\t{}\n'.format(total_jcec_event_ids_no)
        ])
        open(ofile, 'w').writelines(lines)


@pkgsfuncdeco
def export_psi_stats(dct, ofile):
    lines = ['stat_item\tvalue\n']
    total_jc_event_ids_no = 0
    total_jcec_event_ids_no = 0
    s1_jc_in_total = 0
    s1_jc_ex_total = 0
    s1_jc_total = 0
    s1_jcec_in_total = 0
    s1_jcec_ex_total = 0
    s1_jcec_total = 0
    for et, jt2s1dct in dct.items():
        jc_event_ids = set()
        jc_event_ids.update(jt2s1dct['JC']['s1in'], jt2s1dct['JC']['s1ex'])
        jcec_event_ids = set()
        jcec_event_ids.update(jt2s1dct['JCEC']['s1in'], jt2s1dct['JCEC']['s1ex'])
        lines.extend([
            '{}_SAMPLE_1_JunctionCountOnly_inclusion\t{}\n'.format(et, len(jt2s1dct['JC']['s1in'])),
            '{}_SAMPLE_1_JunctionCountOnly_exclusion\t{}\n'.format(et, len(jt2s1dct['JC']['s1ex'])),
            '{}_SAMPLE_1_JunctionCountOnly_total\t{}\n'.format(et, len(jt2s1dct['JC']['s1in']) + len(
                jt2s1dct['JC']['s1ex'])),
            '{}_SAMPLE_1_ReadsOnTargetAndJunctionCounts_inclusion\t{}\n'.format(et, len(jt2s1dct['JCEC']['s1in'])),
            '{}_SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion\t{}\n'.format(et, len(jt2s1dct['JCEC']['s1ex'])),
            '{}_SAMPLE_1_ReadsOnTargetAndJunctionCounts_total\t{}\n'.format(et, len(jt2s1dct['JCEC']['s1in']) + len(
                jt2s1dct['JCEC']['s1ex'])),
            '{}_SAMPLE_2_JunctionCountOnly_inclusion\t{}\n'.format(et, len(jt2s1dct['JC']['s1ex'])),
            '{}_SAMPLE_2_JunctionCountOnly_exclusion\t{}\n'.format(et, len(jt2s1dct['JC']['s1in'])),
            '{}_SAMPLE_2_JunctionCountOnly_total\t{}\n'.format(et, len(jt2s1dct['JC']['s1ex']) + len(
                jt2s1dct['JC']['s1in'])),
            '{}_SAMPLE_2_ReadsOnTargetAndJunctionCounts_inclusion\t{}\n'.format(et, len(jt2s1dct['JCEC']['s1ex'])),
            '{}_SAMPLE_2_ReadsOnTargetAndJunctionCounts_exclusion\t{}\n'.format(et, len(jt2s1dct['JCEC']['s1in'])),
            '{}_SAMPLE_2_ReadsOnTargetAndJunctionCounts_total\t{}\n'.format(et, len(jt2s1dct['JCEC']['s1ex']) + len(
                jt2s1dct['JCEC']['s1in'])),
        ])
        s1_jc_in_total += len(jt2s1dct['JC']['s1in'])
        s1_jc_ex_total += len(jt2s1dct['JC']['s1ex'])
        s1_jc_total = s1_jc_in_total + s1_jc_ex_total
        s1_jcec_in_total += len(jt2s1dct['JCEC']['s1in'])
        s1_jcec_ex_total += len(jt2s1dct['JCEC']['s1ex'])
        s1_jcec_total = s1_jcec_in_total + s1_jcec_ex_total
        total_jc_event_ids_no += len(jc_event_ids)
        total_jcec_event_ids_no += len(jcec_event_ids)
    else:
        lines.extend([
            'SAMPLE_1_JunctionCountOnly_inclusion_total\t{}\n'.format(s1_jc_in_total),
            'SAMPLE_1_JunctionCountOnly_exclusion_total\t{}\n'.format(s1_jc_ex_total),
            'SAMPLE_1_JunctionCountOnly_total\t{}\n'.format(s1_jc_total),
            'SAMPLE_1_ReadsOnTargetAndJunctionCounts_inclusion_total\t{}\n'.format(s1_jcec_in_total),
            'SAMPLE_1_ReadsOnTargetAndJunctionCounts_exclusion_total\t{}\n'.format(s1_jcec_ex_total),
            'SAMPLE_1_ReadsOnTargetAndJunctionCounts_total\t{}\n'.format(s1_jcec_total),
            'SAMPLE_2_JunctionCountOnly_inclusion_total\t{}\n'.format(s1_jc_ex_total),
            'SAMPLE_2_JunctionCountOnly_exclusion_total\t{}\n'.format(s1_jc_in_total),
            'SAMPLE_2_JunctionCountOnly_total\t{}\n'.format(s1_jc_total),
            'SAMPLE_2_ReadsOnTargetAndJunctionCounts_inclusion_total\t{}\n'.format(s1_jcec_ex_total),
            'SAMPLE_2_ReadsOnTargetAndJunctionCounts_exclusion_total\t{}\n'.format(s1_jcec_in_total),
            'SAMPLE_2_ReadsOnTargetAndJunctionCounts_total\t{}\n'.format(s1_jcec_total),
            'total_JunctionCountOnly\t{}\n'.format(total_jc_event_ids_no),
            'total_ReadsOnTargetAndJunctionCounts\t{}\n'.format(total_jcec_event_ids_no),
        ])
        open(ofile, 'w').writelines(lines)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate rmats statistics results')
    parser.add_argument('-i', action='store', required=True,
                        help='downloaded directory', dest='root')
    parser.add_argument('-m', action='store', choices=['pvalue', 'fdr'], required=True,
                        help='testing method', dest='pvalue_fdr')
    parser.add_argument('-s', action='store', required=True, type=str,
                        help='document main id', dest='main_id')
    parser.add_argument('-c', action='store', required=True, type=float,
                        help='significance threshold', dest='fdr')
    parser.add_argument('-p', action='store', required=True, type=float,
                        help='inc level diff threshold', dest='psi')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', dest='output')
    parser.add_argument('-v', action='store', required=True, type=int, default=1,
                        help='mongo_db_version', dest='version')

    args = parser.parse_args()
    main(args)
