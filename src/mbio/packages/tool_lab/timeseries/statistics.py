# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
from optparse import OptionParser

import pandas as pd

parser = OptionParser(description='Conduct statistics of the results of STEM')
parser.add_option('--method', dest='method', choices=['SCM', 'K'], help='clustering method')
parser.add_option('--genetable', dest='genetable', type=str, help='gene table of STEM results')
parser.add_option('-p', dest='p', type=str, help='profile table of STEM results')
parser.add_option('-k', dest='k', type=str, help='kmeans cluster table of STEM results')
parser.add_option('--matrix', dest='matrix', type=str, help='raw data matrix')
parser.add_option('--output', dest='output', type=str, help='output directory')
(opts, args) = parser.parse_args()


def main(opts):
    print 'INFO: start reading {}'.format(opts.genetable)
    detailtable = os.path.join(opts.output, 'detail.{}.tsv'.format(opts.method))
    export_detail_table(opts.method, opts.genetable, opts.matrix, detailtable)
    if os.path.getsize(detailtable):
        print 'INFO: succeed in exporting {}'.format(detailtable)
    classtable = {'SCM': opts.p, 'K': opts.k}[opts.method]
    print 'INFO: start reading {}'.format(classtable)
    clustertable = os.path.join(opts.output, 'cluster.{}.tsv'.format(opts.method))
    export_cluster_table(opts.method, classtable, clustertable)
    if os.path.getsize(clustertable):
        print 'INFO: succeed in exporting {}'.format(clustertable)


def export_detail_table(method, genetable, matrix, detailtable):
    {'SCM': process_scm_detail, 'K': process_k_detail}[method](genetable, matrix, detailtable)


def export_cluster_table(method, genetable, clustertable):
    {'SCM': process_scm_cluster, 'K': process_k_cluster}[method](genetable, clustertable)


def process_scm_detail(genetable, matrix, detailtable):
    lst = list()
    df = pd.read_table(genetable)
    columns = ['seq_id', 'profile'] + list(df.columns)[3:]
    for i, row in df.iterrows():
        for p in str(row['Profile']).split(';'):
            dct = row.to_dict()
            dct.pop('SPOT')
            dct.pop('Profile')
            dct.update({'profile': p})
            lst.append(dct)
    else:
        df = pd.DataFrame(lst, columns=columns)
        df_matrix = pd.read_table(matrix)
        if 'seq_id' in df_matrix.columns.values:
            dct = {s.upper(): s for s in df_matrix['seq_id']}
        elif 'gene_id' in df_matrix.columns.values and len(df_matrix['gene_id']) == len(set(df_matrix['gene_id'])):
            df_matrix.rename(columns={'gene_id': 'seq_id'}, inplace=True)
            dct = {s.upper(): s for s in df_matrix['seq_id']}
        elif 'transcript_id' in df_matrix.columns.values:
            df_matrix.rename(columns={'transcript_id': 'seq_id'}, inplace=True)
            dct = {s.upper(): s for s in df_matrix['seq_id']}
        else:
            raw_id = list(df_matrix.columns.values)[0]
            df_matrix.rename(columns={raw_id: 'seq_id'}, inplace=True)
            dct = {s.upper(): s for s in df_matrix['seq_id']}
            # raise Exception("Can't find seq_id!")
        for i in df.index:
            df.loc[i, 'seq_id'] = dct[df.loc[i, 'seq_id']]
        else:
            df.to_csv(detailtable, sep='\t', index=False)


def process_k_detail(genetable, matrix, detailtable):
    lst = list()
    df = pd.read_table(genetable)
    df = df.rename(columns={'Cluster': 'cluster'})
    columns = ['seq_id', 'cluster'] + list(df.columns)[3:]
    df = df.reindex(columns, axis=1)
    df_matrix = pd.read_table(matrix)
    if 'seq_id' in df_matrix.columns.values:
        dct = {s.upper(): s for s in pd.read_table(matrix)['seq_id']}
    elif 'gene_id' in df_matrix.columns.values and len(df_matrix['gene_id']) == len(set(df_matrix['gene_id'])):
        df_matrix.rename(columns={'gene_id': 'seq_id'}, inplace=True)
        dct = {s.upper(): s for s in df_matrix['seq_id']}
    elif 'transcript_id' in df_matrix.columns.values:
        df_matrix.rename(columns={'transcript_id': 'seq_id'}, inplace=True)
        dct = {s.upper(): s for s in df_matrix['seq_id']}
    for i in df.index:
        df.loc[i, 'seq_id'] = dct[df.loc[i, 'seq_id']]
    else:
        df.to_csv(detailtable, sep='\t', index=False)


def process_scm_cluster(ptable, clustertable):
    df = pd.read_table(ptable)
    columns = {'Profile ID': 'profile', 'Profile Model': 'model', 'Cluster (-1 non-significant)': 'cluster',
               '# Genes Assigned': 'assigned', '# Gene Expected': 'expected', 'p-value': 'pvalue'}
    df = df.rename(columns=columns)
    df.to_csv(clustertable, sep='\t', index=False)


def process_k_cluster(ktable, clustertable):
    df = pd.read_table(ktable)
    columns = {'Cluster': 'cluster', 'Cluster Mean': 'model', 'Number of Genes': 'number'}
    df = df.rename(columns=columns)
    df.to_csv(clustertable, sep='\t', index=False)


if __name__ == '__main__':
    if all(map(hasattr, [opts] * 4, ['method', 'genetable', 'output', 'matrix'])):
        main(opts)
    else:
        parser.print_help()
