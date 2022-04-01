import pandas as pd
import dask
import glob
import os
import re
import numpy as np
import json

def volcano(detail_df_dict, diff_output, diff_method, output_dir, pvalue_padjust=None):
    volcano_list = list()
    sig_status = dict()
    volcano_df_dict = dict()
    for (ctrl, test), diff_pd in detail_df_dict.items():
        status_list = _get_volcano_status_cutoff(diff_pd, pvalue_padjust)
        cmp_combine = ctrl + '|' + test
        sig_status[cmp_combine] = status_list
        cmp_pd = pd.DataFrame([cmp_combine] * diff_pd.shape[0], columns=['compare'])
        if diff_method.lower() in ['deseq2', 'edger', 'degseq', 'limma', 'svaseqlimma']:
            volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', pvalue_padjust, 'significant', 'regulate']]
            bool_ind = volcano_pd[pvalue_padjust] <= 0
            min_pvalue = min([x if x > 0 else '' for x in volcano_pd[pvalue_padjust].tolist()])
            volcano_pd.loc[bool_ind, pvalue_padjust] = min_pvalue
            volcano_pd[pvalue_padjust] = -np.log10(volcano_pd[pvalue_padjust])
            volcano_pd.dropna(inplace=True)
            volcano_pd.columns = ['seq_id', 'log2fc', 'log10pvalue', 'significant', 'regulate']
            # bool_ind = volcano_pd['log10pvalue'] > stat_cutoff
        if diff_method.lower() in ['noiseq']:
            volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', 'D', 'significant', 'regulate']]
            volcano_pd.dropna(inplace=True)
            volcano_pd.columns = ['seq_id', 'log2fc', 'D', 'significant', 'regulate']
        volcano_pd = pd.concat([volcano_pd, cmp_pd], axis=1)
        volcano_pd_nosig = volcano_pd[volcano_pd['significant'] == 'no']
        # random select 10000 not sig diff genes for plotting
        if volcano_pd_nosig.shape[0] > 10000:
            volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
            volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.5)
            while volcano_pd_nosig.shape[0] > 50000:
                volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
            volcano_pd = pd.concat([volcano_pd_sig, volcano_pd_nosig], axis=0)
        volcano_list.append(volcano_pd)
        volcano_df_dict[(ctrl, test)] = volcano_pd
    all_volcano = reduce(lambda x, y: pd.concat([x, y]), volcano_list)
    all_volcano.to_csv(os.path.join(output_dir, 'all_volcano.txt'), sep='\t', index=False)
    return volcano_df_dict, sig_status

def scatter(detail_tmp_dict, volcano_df_dict, diff_output, diff_method, output_dir):
    scatter_list = list()
    scatter_df_dict = dict()
    for (ctrl, test), tmp_pd in detail_tmp_dict.items():
        volcano_pd = volcano_df_dict[(ctrl, test)]
        if diff_method.lower() in ['deseq2', 'edger', 'degseq', 'limma', 'svaseqlimma']:
            scatter_pd = tmp_pd.loc[:, ['seq_id', 'group1', 'group2', 'compare', 'significant', 'regulate']]
            scatter_pd.set_index('seq_id', inplace=True)
            scatter_pd = scatter_pd.loc[volcano_pd['seq_id'], :].reset_index()
            scatter_pd['group1'] = (scatter_pd['group1'] + 1).apply(np.log10)
            scatter_pd['group2'] = (scatter_pd['group2'] + 1).apply(np.log10)
        if diff_method.lower() in ['noiseq']:
            scatter_pd = tmp_pd.loc[:, ['seq_id', 'compare', 'significant', 'regulate', '{}_mean'.format(ctrl), '{}_mean'.format(test)]]
            scatter_pd.set_index('seq_id', inplace=True)
            scatter_pd = scatter_pd.loc[volcano_pd['seq_id'], :].reset_index()
            scatter_pd['group1'] = (scatter_pd['{}_mean'.format(ctrl)]+1).apply(np.log10)
            scatter_pd['group2'] = (scatter_pd['{}_mean'.format(test)]+1).apply(np.log10)
        scatter_list.append(scatter_pd)
        scatter_df_dict[(ctrl, test)] = scatter_pd
    all_scatter = reduce(lambda x, y: pd.concat([x, y]), scatter_list)
    all_scatter.to_csv(os.path.join(output_dir, 'all_scatter.txt'), sep='\t', index=False)

def detail(diff_output, diff_method, output_dir):
    diff_files = glob.glob(os.path.join(diff_output, '*_vs_*.*.xls'))
    diff_list = list()
    cmp_list = list()
    cmp_detail_dict = dict()
    detail_df_dict = dict()
    detail_tmp_dict = dict()
    final_details = []
    for each in diff_files:
        if each.endswith('.annot.xls'):
            continue
        if each.endswith('.normalize.xls'):
            continue
        if each.endswith('.sizeFactor.xls'):
            continue
        if each.endswith('.normFactor.xls'):
            continue
        final_detail_df = dask.delayed(detail_extract)(each, diff_method)
        final_details.append(final_detail_df)
    else:
        final_details = dask.compute(*final_details)
        for ctrl, test, samples, tmp_pd, diff_pd in final_details:
            diff_list.append(tmp_pd)
            detail_df_dict[(ctrl, test)] = diff_pd
            detail_tmp_dict[(ctrl, test)] = tmp_pd
            cmp_combine = ctrl + '|' + test
            cmp_list.append(cmp_combine)
            cmp_detail_dict[cmp_combine] = samples
        all_detail = reduce(lambda x, y: pd.concat([x, y]), diff_list)
        all_detail.to_csv(os.path.join(output_dir, 'all_detail.txt'), sep='\t', index=False)
        return detail_df_dict, detail_tmp_dict, cmp_list, cmp_detail_dict
    #
    #     diff_pd = pd.read_table(each, header=0, sep='\t')
    #     columns = diff_pd.columns
    #     fname = os.path.basename(each)
    #     ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(diff_method.lower()), fname).groups()
    #     fc_ind = list(columns).index('fc')
    #     if diff_method.lower() in ['deseq2', 'edger', 'degseq', 'limma']:
    #         need_cols = ['seq_id', 'fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
    #         need_cols += [columns[fc_ind - 2], columns[fc_ind - 1]]
    #     if diff_method.lower() in ['noiseq']:
    #         need_cols = ['seq_id', '{}_mean'.format(ctrl), '{}_mean'.format(test), 'fc', 'log2fc', 'D', 'prob',
    #                      'significant', 'regulate']
    #         need_cols += [columns[fc_ind - 4], columns[fc_ind - 3]]
    #     print(need_cols)
    #     samples = list()
    #     for x in columns:
    #         _m = re.match(r'(.*)_count$', x)
    #         if _m:
    #             samples.append(_m.groups()[0])
    #
    #     cmp_combine = ctrl + '|' + test
    #     cmp_list.append(cmp_combine)
    #     cmp_detail_dict[cmp_combine] = samples
    #     cmp_pd = pd.DataFrame([cmp_combine] * diff_pd.shape[0], columns=['compare'])
    #     tmp_pd = pd.concat([diff_pd.loc[:, need_cols], cmp_pd], axis=1)
    #     tmp_pd.columns = list(tmp_pd.columns[:-3]) + ['group1', 'group2', 'compare']
    #     print(tmp_pd.columns)
    #     print(tmp_pd)
    #     diff_list.append(tmp_pd)
    #     detail_df_dict[(ctrl, test)] = diff_pd
    #     detail_tmp_dict[(ctrl, test)] = tmp_pd
    # all_detail = reduce(lambda x, y: pd.concat([x, y]), diff_list)
    # all_detail.to_csv(os.path.join(output_dir, 'all_detail.txt'), sep='\t', index=False)
    # return detail_df_dict, detail_tmp_dict, cmp_list, cmp_detail_dict

def detail_extract(each, diff_method):
    diff_pd = pd.read_table(each, header=0, sep='\t')
    columns = diff_pd.columns
    fname = os.path.basename(each)
    ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(diff_method.lower()), fname).groups()
    fc_ind = list(columns).index('fc')
    if diff_method.lower() in ['deseq2', 'edger', 'degseq', 'limma', 'svaseqlimma']:
        need_cols = ['seq_id', 'fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
        need_cols += [columns[fc_ind - 2], columns[fc_ind - 1]]
    if diff_method.lower() in ['noiseq']:
        need_cols = ['seq_id', '{}_mean'.format(ctrl), '{}_mean'.format(test), 'fc', 'log2fc', 'D', 'prob',
                     'significant', 'regulate']
        need_cols += [columns[fc_ind - 4], columns[fc_ind - 3]]
    samples = list()
    for x in columns:
        _m = re.match(r'(.*)_count$', x)
        if _m:
            samples.append(_m.groups()[0])
    cmp_combine = ctrl + '|' + test
    cmp_pd = pd.DataFrame([cmp_combine] * diff_pd.shape[0], columns=['compare'])
    tmp_pd = pd.concat([diff_pd.loc[:, need_cols], cmp_pd], axis=1)
    tmp_pd.columns = list(tmp_pd.columns[:-3]) + ['group1', 'group2', 'compare']
    return ctrl, test, samples, tmp_pd, diff_pd

def _get_volcano_status_cutoff(diff_table, pvalue_padjust):
    sig_status = list()
    sig_mark = diff_table['significant']
    reg_list = diff_table['regulate']
    if 'no' in list(sig_mark):
        no_sig_num = sig_mark[sig_mark == "no"].shape[0]
        sig_status.append('nosig_'+str(no_sig_num))
    if 'yes' in list(sig_mark):
        reg_mark = reg_list[sig_mark == 'yes']
        if 'down' in list(reg_mark):
            down_num = reg_mark[reg_mark == 'down'].shape[0]
            sig_status.append('down_'+str(down_num))
        if 'up' in list(reg_mark):
            up_num = reg_mark[reg_mark == 'up'].shape[0]
            sig_status.append('up_'+str(up_num))

    # sig_pvalues = diff_table[pvalue_padjust][diff_table['significant'] == "yes"]
    # log10_sig_pvalues = -np.log10(sig_pvalues)
    # log10_pvalue_list = sorted(list(log10_sig_pvalues[log10_sig_pvalues > 0]))
    #
    # if len(sig_pvalues) > 2000:
    #     log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.85)]
    # elif len(sig_pvalues) > 1000:
    #     log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.90)]
    # elif len(sig_pvalues) > 500:
    #     log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.95)]
    # elif len(sig_pvalues) > 250:
    #     log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.99)]
    # elif len(sig_pvalues) == 0:
    #     tmp = -np.log10(diff_table[pvalue_padjust])
    #     tmp_list = sorted(tmp[tmp > 0])
    #     if len(tmp_list) == 0:
    #         log10_pvalue_cutoff = 200
    #     else:
    #         log10_pvalue_cutoff = tmp_list[int(len(tmp_list)*0.9)]
    # else:
    #     # print(pvalue_padjust, diff_table, log10_pvalue_list)
    #     log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list)*0.8)]
    return sig_status

def write_json(cmp_list, cmp_detail_dict, sig_status, output_dir):
    json.dump({
        'cmp_list': cmp_list,
        'cmp_detail_dict': cmp_detail_dict,
        'sig_status': sig_status
    }, open(os.path.join(output_dir, 'json'), 'w'), indent=4)


def main(args):
    diff_method = args.program
    input_dir = args.input_dir
    output_dir = args.output_dir
    detail_df_dict, detail_tmp_dict,  cmp_list, cmp_detail_dict = detail(input_dir, diff_method, output_dir)
    if diff_method.lower() in ['deseq2', 'degseq', 'edger', 'limma', 'svaseqlimma']:
        volcano_df_dict, sig_status = volcano(detail_df_dict, input_dir, diff_method, output_dir, pvalue_padjust=args.pvalue_padjust)
    if diff_method.lower() in ['noiseq']:
        volcano_df_dict, sig_status = volcano(detail_df_dict, input_dir, diff_method, output_dir)
    scatter(detail_tmp_dict, volcano_df_dict, input_dir, diff_method, output_dir)
    write_json(cmp_list, cmp_detail_dict, sig_status, output_dir)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate files for dumping data into mongo')

    parser.add_argument('-m', action='store', choices=['DESeq2', 'edgeR', 'DEGseq', 'Limma', 'NOIseq', 'svaseqlimma'], required=True,
                        help='selected program', dest='program')
    parser.add_argument('-i', action='store', required=True,
                        help='input directory', metavar='<DIR>', dest='input_dir')
    parser.add_argument('-o', action='store', required=True,
                        help='input directory', metavar='<DIR>', dest='output_dir')
    parser.add_argument('-p', action='store', choices=['pvalue', 'padjust'], required=True,
                        help='statistics type', dest='pvalue_padjust')

    args = parser.parse_args()

    main(args)
