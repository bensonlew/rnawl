# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import pickle
import subprocess

import numpy as np
import pandas as pd


def main(args):
    df = pd.read_table(args.matrix)
    seq_ids = list({line.strip().split('\t')[0] for n, line in enumerate(open(args.result)) if n != 0})
    desdf = pd.read_table(args.design)
    for group in desdf.columns[3:]:
        subdesdf = desdf[desdf[group] == 1]
        subdesdf = subdesdf.sort_values(['Replicate', 'Sample'])
        subdesdf = subdesdf.reindex(['Sample', 'Time', 'Replicate', group], axis=1)
        design = os.path.join(os.path.dirname(args.output), 'design.{}.txt'.format(group))
        subdesdf.to_csv(design, sep='\t', index=False)
        samples = list(subdesdf['Sample'])
        subdf = df.rename({{'T': 'transcript_id', 'G': 'gene_id'}[args.level]: 'seq_id'}, axis=1)
        subdf = subdf.reindex(['seq_id'] + samples, axis=1)
        subdf = subdf.query('seq_id in @seq_ids')
        subdf, data = time_sort(subdf, subdesdf)
        pickle.dump(data, open(os.path.join(args.output, '{}.pkl'.format(group)), 'w'))
        heatmap = os.path.join(os.path.dirname(args.output), 'heatmap.{}.tsv'.format(group))
        subdf.to_csv(heatmap, sep='\t')
        cmd = '{} {} -i {} -d {} -n {} -t {} -o {}'.format(
            args.interpreter, args.script, heatmap, design, group, args.type.upper(), args.output
        )
        proc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        retcode = proc.wait()
        outs, errs = proc.communicate()
        print 'INFO: stdout of subprocess ->\n{}'.format(outs)
        print 'INFO: stderr of subprocess ->\n{}'.format(errs)
        if retcode:
            raise Exception('ERROR: fail to excecute command -> {} (returncode {})'.format(cmd, retcode))
        else:
            print('INFO: succeed in executing command -> {}'.format(cmd))


def time_sort(raw_df, des_df, timeseries=None):
    num_df = raw_df.set_index('seq_id')
    scale_row = lambda row: (row - np.min(row)) / (np.max(row) - np.min(row))
    scale_df = pd.DataFrame(scale_row(row) for i, row in num_df.iterrows())
    index = scale_df.index
    des_df = des_df.sort_values('Time')
    groups = des_df['Time'].tolist()
    samples = des_df['Time'].tolist()
    columns = pd.MultiIndex.from_arrays([groups, samples], names=['group', 'sample'])
    df = pd.DataFrame(np.array(scale_df), index=index, columns=columns)
    gmdf = df.groupby(level='group', axis=1).mean()
    dct = {c: list() for c in gmdf.columns}
    for i, row in gmdf.iterrows():
        if row.idxmax() in dct:
            dct[row.idxmax()].append(row)
    c2i = {c: pd.DataFrame(dct[c]).sort_values(c, ascending=False).index for c in gmdf.columns if dct[c]}
    ret_df = pd.concat(scale_df.reindex(c2i[c]) for c in gmdf.columns if c in c2i)
    ret_df.index.name = 'seq_id'
    ret_df = ret_df.reset_index()
    ret_df = ret_df.drop_duplicates('seq_id').set_index('seq_id')
    data = list()
    for time, rows in dct.items():
        seq_ids = list({row.name for row in rows})
        document = {'time': time, 'seq_ids': seq_ids}
        data.append(document)
    return ret_df, data


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Obtain data for drawing heatmap')
    parser.add_argument('-m', dest='matrix', type=str, help='input raw data matrix table', required=True)
    parser.add_argument('-d', dest='design', type=str, help='input design table', required=True)
    parser.add_argument('-r', dest='result', type=str, help='input masigpro result table', required=True)
    parser.add_argument('--interpreter', dest='interpreter', type=str, help='path of Rscript interpreter',
                        required=True)
    parser.add_argument('--script', dest='script', type=str, help='path of script for drawing heatmap', required=True)
    parser.add_argument('-l', dest='level', type=str, help='expression level', required=True)
    parser.add_argument('-t', dest='type', type=str, help='expression way', required=True)
    parser.add_argument('-o', dest='output', type=str, help='output directory for heatmap', required=True)
    args = parser.parse_args()

    main(args)
