# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os

import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


def main(args):
    program = args.program
    input_dir = args.input_dir
    exp_matrix, group_table, control_table, kind_table = \
        args.exp_matrix, args.group_table, args.control_table, args.kind_table
    filter, threshold = args.filter, args.threshold
    method = args.method
    stat_type, stat_cutoff, fc = args.stat_type, args.stat_cutoff, args.fc
    output_dir = args.output_dir

    exp_df, group_dict, control_list, kind_df = parse_input(exp_matrix, group_table, control_table, kind_table)
    logging.info('succeed in parsing input files')
    remain_ids_dict = retain_ids(exp_df, group_dict, control_list, filter, threshold)
    logging.info('succeed in collecting ids by screening')
    if program == 'DESeq2':
        detail_df_dict = call_detail_deseq2(input_dir, exp_df, group_dict, control_list, remain_ids_dict, method,
                                            stat_type, stat_cutoff, fc, kind_df, output_dir)
    elif program == 'edgeR':
        detail_df_dict = call_detail_edger(input_dir, exp_df, group_dict, control_list, remain_ids_dict, method,
                                           stat_type, stat_cutoff, fc, kind_df, output_dir)
    elif program == 'DEGseq':
        detail_df_dict = call_detail_degseq(input_dir, exp_df, group_dict, control_list, remain_ids_dict, method,
                                            stat_type, stat_cutoff, fc, kind_df, output_dir)
    logging.info('succeed in exporting detail data from {} result'.format(program))
    call_summary(output_dir, control_list, remain_ids_dict, detail_df_dict)
    logging.info('succeed in exporting summary data')
    call_volcano(detail_df_dict, stat_type, output_dir)
    logging.info('succeed in exporting volcano data')
    call_scatter(detail_df_dict, exp_df, group_dict, output_dir)
    logging.info('succeed in exporting scatter data')


def call_scatter(detail_df_dict, exp_df, group_dict, output_dir):
    scatter_df_dict = dict()
    for (ctrl, case), detail_df in detail_df_dict.items():
        scatter_df = detail_df.reindex(['compare', 'significant', 'regulate'], axis=1)
        scatter_df['group1'] = exp_df.reindex(group_dict[ctrl], axis=1).mean(axis=1).apply(lambda x: np.log10(x + 1))
        scatter_df['group2'] = exp_df.reindex(group_dict[case], axis=1).mean(axis=1).apply(lambda x: np.log10(x + 1))
        scatter_df_dict[(ctrl, case)] = scatter_df
    else:
        all_scatter_df = pd.concat(scatter_df_dict.values())
        all_scatter_df.to_csv(os.path.join(output_dir, 'scatter.txt'), sep='\t')


def call_volcano(detail_df_dict, stat_type, output_dir):
    volcano_df_dict = dict()
    for (ctrl, case), detail_df in detail_df_dict.items():
        detail_df['log10pvalue'] = -np.log10(detail_df[stat_type])
        volcano_df = detail_df.reset_index()
        volcano_df = volcano_df.reindex(['seq_id', 'log2fc', 'log10pvalue', 'significant', 'regulate', 'compare'],
                                        axis=1)
        # issig_volcano_df = volcano_df[volcano_df['significant'] == 'yes']
        # nosig_volcano_df = volcano_df[volcano_df['significant'] == 'no']
        # if nosig_volcano_df.shape[0] > 5000:
        #     frac = 5000.0 / nosig_volcano_df.shape[0]
        #     nosig_volcano_df = nosig_volcano_df.sample(frac=frac)
        # volcano_df = pd.concat([issig_volcano_df, nosig_volcano_df])
        volcano_df_dict[(ctrl, case)] = volcano_df

    all_volcano_df = pd.concat(volcano_df_dict.values())
    max_log10value = max(all_volcano_df[all_volcano_df['log10pvalue'] != np.inf]['log10pvalue'])
    all_volcano_df['log10pvalue'] = all_volcano_df['log10pvalue'].apply(lambda x: max_log10value if x == np.inf else x)
    all_volcano_df.to_csv(os.path.join(output_dir, 'volcano.txt'), sep='\t', index=False)


def call_summary(output_dir, control_list, remain_ids_dict, detail_df_dict):
    remain_id_set = set()
    remain_id_set.update(*remain_ids_dict.values())
    sig_id_dict = dict()
    for (ctrl, case), detail_df in detail_df_dict.items():
        sig_id_set = set(detail_df[detail_df['significant'] == 'yes'].index)
        sig_id_dict['{}_vs_{}'.format(ctrl, case)] = sig_id_set
    summary_data = list()
    for seq_id in remain_id_set:
        summary_dict = {'seq_id': seq_id, 'sum': 0}
        for pair, sig_id_set in sig_id_dict.items():
            ctrl, case = pair.split('_vs_')
            if seq_id in sig_id_set:
                summary_dict[pair] = 'yes|{}'.format(detail_df_dict[(ctrl, case)].loc[seq_id, 'regulate'])
                summary_dict['sum'] += 1
            else:
                if seq_id in detail_df_dict[(ctrl, case)].index:
                    summary_dict[pair] = 'no|{}'.format(detail_df_dict[(ctrl, case)].loc[seq_id, 'regulate'])
                else:
                    summary_dict[pair] = 'no|no_test'
        if summary_dict['sum']:
            summary_data.append(summary_dict)
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.reindex(['seq_id'] + sig_id_dict.keys() + ['sum'], axis=1)
    summary_df.to_csv(os.path.join(output_dir, 'summary.txt'), sep='\t', index=False)


def call_detail_deseq2(input_dir, exp_df, group_dict, control_list, remain_ids_dict, method, stat_type, stat_cutoff, fc,
                       kind_df, output_dir):
    detail_df_dict = dict()
    for pair_dict in control_list:
        ctrl, case = pair_dict['ctrl'], pair_dict['case']
        remain_ids = remain_ids_dict[(ctrl, case)]
        deseq2_table = os.path.join(input_dir, '{}_vs_{}.deseq2.txt'.format(ctrl, case))
        deseq2_df = pd.read_table(deseq2_table, index_col=0)
        deseq2_df.index.name = 'seq_id'
        deseq2_df = deseq2_df.reindex(remain_ids)
        deseq2_df['log2fc'] = deseq2_df['log2FoldChange'].fillna(0)
        deseq2_df['fc'] = 2 ** deseq2_df['log2fc']
        if method.lower() == 'bh':
            deseq2_df.rename({'padj': 'padjust'}, axis=1, inplace=True)
        else:
            deseq2_df['padjust'] = pd.Series(run_multitest(deseq2_df['pvalue'], method), index=deseq2_df.index)
        deseq2_df['significant'] = ['yes' if i and j else 'no' for i, j in
                                    zip((deseq2_df['fc'].fillna(1) >= fc) | (deseq2_df['fc'].fillna(1) <= 1 / fc),
                                        deseq2_df[stat_type].fillna(1) <= stat_cutoff)]
        deseq2_df['regulate'] = deseq2_df['log2fc'].fillna(0).apply(
            lambda x: 'up' if x > 0 else 'down' if x else 'no_change')
        deseq2_df['compare'] = '{}|{}'.format(ctrl, case)
        detail_df = deseq2_df.reindex(['fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate', 'compare'],
                                      axis=1)
        detail_df['group1'] = exp_df.reindex(group_dict[ctrl], axis=1).mean(axis=1)
        detail_df['group2'] = exp_df.reindex(group_dict[case], axis=1).mean(axis=1)
        na_index = deseq2_df[pd.isna(deseq2_df['log2FoldChange'])].index
        detail_df['pvalue'] = detail_df['pvalue'].fillna(1)
        detail_df['padjust'] = detail_df['padjust'].fillna(1)
        detail_df.loc[na_index, 'significant'] = 'no_test'
        detail_df.loc[na_index, 'regulate'] = 'no_test'
        detail_df = detail_df.join(kind_df)
        detail_df.to_csv(os.path.join(output_dir, '{}_vs_{}.detail.txt'.format(ctrl, case)), sep='\t')
        detail_df_dict[(ctrl, case)] = detail_df
    else:
        return detail_df_dict


def call_detail_edger(input_dir, exp_df, group_dict, control_list, remain_ids_dict, method, stat_type, stat_cutoff, fc,
                      kind_df, output_dir):
    detail_df_dict = dict()
    for pair_dict in control_list:
        ctrl, case = pair_dict['ctrl'], pair_dict['case']
        remain_ids = remain_ids_dict[(ctrl, case)]
        edger_table = os.path.join(input_dir, '{}_vs_{}.edger.txt'.format(ctrl, case))
        edger_df = pd.read_table(edger_table, index_col=0)
        edger_df.index.name = 'seq_id'
        edger_df = edger_df.reindex(remain_ids)
        edger_df['log2fc'] = edger_df['logFC'].fillna(0)
        edger_df['fc'] = 2 ** edger_df['log2fc']
        edger_df['pvalue'] = edger_df['PValue']
        edger_df['padjust'] = pd.Series(run_multitest(edger_df['pvalue'], method), index=edger_df.index)
        edger_df['significant'] = ['yes' if i and j else 'no' for i, j in
                                   zip((edger_df['fc'].fillna(1) >= fc) | (edger_df['fc'].fillna(1) <= 1 / fc),
                                       edger_df[stat_type].fillna(1) <= stat_cutoff)]
        edger_df['regulate'] = edger_df['log2fc'].fillna(0).apply(
            lambda x: 'up' if x > 0 else 'down' if x else 'no_change')
        edger_df['compare'] = '{}|{}'.format(ctrl, case)
        detail_df = edger_df.reindex(['fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate', 'compare'],
                                     axis=1)
        detail_df['group1'] = exp_df.reindex(group_dict[ctrl], axis=1).mean(axis=1)
        detail_df['group2'] = exp_df.reindex(group_dict[case], axis=1).mean(axis=1)
        na_index = edger_df[pd.isna(edger_df['logFC'])].index
        detail_df['pvalue'] = detail_df['pvalue'].fillna(1)
        detail_df['padjust'] = detail_df['padjust'].fillna(1)
        detail_df.loc[na_index, 'significant'] = 'no_test'
        detail_df.loc[na_index, 'regulate'] = 'no_test'
        detail_df = detail_df.join(kind_df)
        detail_df.to_csv(os.path.join(output_dir, '{}_vs_{}.detail.txt'.format(ctrl, case)), sep='\t')
        detail_df_dict[(ctrl, case)] = detail_df
    else:
        return detail_df_dict


def call_detail_degseq(input_dir, exp_df, group_dict, control_list, remain_ids_dict, method, stat_type, stat_cutoff, fc,
                       kind_df, output_dir):
    detail_df_dict = dict()
    for pair_dict in control_list:
        ctrl, case = pair_dict['ctrl'], pair_dict['case']
        remain_ids = remain_ids_dict[(ctrl, case)]
        degseq_table = os.path.join(input_dir, '{}_vs_{}/output_score.txt'.format(ctrl, case))
        degseq_df = pd.read_table(degseq_table, index_col=0)
        degseq_df.index.name = 'seq_id'
        degseq_df = degseq_df.reindex(remain_ids)
        degseq_df['log2fc'] = degseq_df['log2(Fold_change) normalized'].fillna(0)
        degseq_df['fc'] = 2 ** degseq_df['log2fc']
        degseq_df['pvalue'] = degseq_df['p-value']
        degseq_df['padjust'] = pd.Series(run_multitest(degseq_df['pvalue'], method), index=degseq_df.index)
        degseq_df['significant'] = ['yes' if i and j else 'no' for i, j in
                                    zip((degseq_df['fc'].fillna(1) >= fc) | (degseq_df['fc'].fillna(1) <= 1 / fc),
                                        degseq_df[stat_type].fillna(1) <= stat_cutoff)]
        degseq_df['regulate'] = degseq_df['log2fc'].fillna(0).apply(
            lambda x: 'up' if x > 0 else 'down' if x else 'no_change')
        degseq_df['compare'] = '{}|{}'.format(ctrl, case)
        detail_df = degseq_df.reindex(['fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate', 'compare'],
                                      axis=1)
        detail_df['group1'] = exp_df.reindex(group_dict[ctrl], axis=1).mean(axis=1)
        detail_df['group2'] = exp_df.reindex(group_dict[case], axis=1).mean(axis=1)
        na_index = degseq_df[pd.isna(degseq_df['log2(Fold_change) normalized'])].index
        detail_df['pvalue'] = detail_df['pvalue'].fillna(1)
        detail_df['padjust'] = detail_df['padjust'].fillna(1)
        detail_df.loc[na_index, 'significant'] = 'no_test'
        detail_df.loc[na_index, 'regulate'] = 'no_test'
        detail_df = detail_df.join(kind_df)
        detail_df.index.name = 'seq_id'
        detail_df.to_csv(os.path.join(output_dir, '{}_vs_{}.detail.txt'.format(ctrl, case)), sep='\t')
        detail_df_dict[(ctrl, case)] = detail_df
    else:
        return detail_df_dict


def run_multitest(pvals, method):
    method = method.lower()
    assert method in ['bonferroni', 'holm', 'bh', 'by']
    method = 'fdr_{}'.format(method) if method in ['bh', 'by'] else method
    pvals = pd.Series(np.array(pvals))
    pvals_notna = pvals[pd.notna(pvals)]
    pvals_isna = pvals[pd.isna(pvals)]
    # reject_, pvals_corrected_, alphacSidak, alphacBonf = multipletests(pvals, method=method)
    pvals_corrected = multipletests(pvals_notna, method=method)[1]
    padjs_notna = pd.Series(pvals_corrected, index=pvals_notna.index)
    padjs_isna = pd.Series(np.nan, index=pvals_isna.index)
    padjs = pd.concat([padjs_notna, padjs_isna])
    padjs = padjs.sort_index()
    return np.array(padjs)


def parse_input(exp_matrix, group_table, control_table, kind_table):
    exp_df = pd.read_table(exp_matrix, index_col=0)
    exp_df.index.name = 'seq_id'
    exp_df = exp_df.reset_index().drop_duplicates().set_index('seq_id')
    group_dict = dict()
    sample_list = list()
    for line in open(group_table):
        if line.strip() and line[0] != '#':
            sample, group = line.strip().split('\t')
            sample_list.append(sample)
            if group in group_dict:
                group_dict[group].append(sample)
            else:
                group_dict[group] = [sample]
    exp_df = exp_df.reindex(sample_list, axis=1)
    control_list = list()
    for line in open(control_table):
        if line.strip() and line[0] != '#':
            ctrl, case = line.strip().split('\t')
            assert ctrl in group_dict
            assert case in group_dict
            control_list.append({'ctrl': ctrl, 'case': case})
    kind_df = pd.read_table(kind_table, index_col=0)
    kind_df.index.name = 'seq_id'
    return exp_df, group_dict, control_list, kind_df


def process_exp_df(exp_df, group_dict=None, log10=True):
    if group_dict:
        group_exp_dfs = list()
        for group in group_dict:
            group_exp_df = exp_df.reindex(group_dict[group], axis=1).mean(axis=1)
            group_exp_df.name = group
            group_exp_dfs.append(group_exp_df)
        else:
            exp_df = pd.concat(group_exp_dfs, axis=1)
    if log10:
        exp_df = np.log10(exp_df)
    return exp_df


def retain_ids(exp_df, group_dict, control_list, filter, threshold):
    remain_ids_dict = dict()
    for pair_dict in control_list:
        ctrl, case = pair_dict['ctrl'], pair_dict['case']
        sample_list = group_dict[ctrl] + group_dict[case]
        use_exp_df = exp_df.reindex(sample_list, axis=1)
        if filter == 'none':
            bool_series = use_exp_df.index == use_exp_df.index
        elif filter == 'mean':
            bool_series = use_exp_df.sum(axis=1) / len(sample_list) > threshold
        elif filter == 'max':
            bool_series = use_exp_df.max(axis=1) > threshold
        elif filter == 'min':
            bool_series = use_exp_df.min(axis=1) > threshold
        elif filter == 'sum':
            bool_series = use_exp_df.sum(axis=1) > threshold
        elif filter == 'median':
            bool_series = use_exp_df.median(axis=1) > threshold
        remain_ids = list(use_exp_df[bool_series].index.drop_duplicates())
        remain_ids_dict[(ctrl, case)] = remain_ids
    else:
        return remain_ids_dict


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate files for dumping data into mongo')

    parser.add_argument('-u', action='store', choices=['DESeq2', 'edgeR', 'DEGseq'], required=True,
                        help='selected program', dest='program')
    parser.add_argument('-i', action='store', required=True,
                        help='input directory', metavar='<DIR>', dest='input_dir')
    parser.add_argument('-e', action='store', required=True,
                        help='expression matrix', metavar='<FILE>', dest='exp_matrix')
    parser.add_argument('-g', action='store', required=True,
                        help='group table', metavar='<FILE>', dest='group_table')
    parser.add_argument('-c', action='store', required=True,
                        help='control table', metavar='<FILE>', dest='control_table')
    parser.add_argument('-k', action='store', required=True,
                        help='kind table', metavar='<FILE>', dest='kind_table')
    parser.add_argument('-f', action='store', choices=['none', 'mean', 'max', 'min', 'sum', 'median'], required=True,
                        help='filter methed', dest='filter')
    parser.add_argument('-t', action='store', type=float, required=True,
                        help='threshold value', metavar='<FLOAT>', dest='threshold')
    parser.add_argument('-m', action='store', choices=['bonferroni', 'holm', 'bh', 'by'], required=True,
                        help='multitest methed', dest='method')
    parser.add_argument('-s', action='store', choices=['pvalue', 'padjust'], required=True,
                        help='statistics type', dest='stat_type')
    parser.add_argument('-p', action='store', type=float, required=True,
                        help='statistics cutoff', metavar='<FLOAT>', dest='stat_cutoff')
    parser.add_argument('-d', action='store', type=float, required=True,
                        help='fold change', metavar='<FLOAT>', dest='fc')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', metavar='<DIR>', dest='output_dir')

    args = parser.parse_args()

    main(args)
