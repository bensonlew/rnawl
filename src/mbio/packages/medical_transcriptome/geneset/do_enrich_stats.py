# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.packages.ref_rna_v2.functions import pkgsfuncdeco
import pandas as pd
import math

@pkgsfuncdeco
def main(args):
    df = pd.read_table(args.ifile)
    df = df[df['enrichment'] != 'p']
    df = df.rename(columns={
        'GO': 'go_id', 'NS': 'go_type', 'name': 'discription', df.columns[9]: 'p_corrected', 'diff_gene': 'seq_list'
    })
    df['study_count'] = df['ratio_in_study'].apply(lambda x: float(x.split('/')[0]))
    df['pop_count'] = df['ratio_in_pop'].apply(lambda x: float(x.split('/')[0]))
    df['enrich_factor'] = df['study_count'] / df['pop_count']
    df['neg_log10p_uncorrected'] = - df['p_uncorrected'].apply(math.log10)
    df['neg_log10p_corrected'] = - df['p_corrected'].apply(math.log10)
    columns = ['go_id', 'go_type', 'enrichment', 'discription', 'ratio_in_study', 'ratio_in_pop', 'p_uncorrected',
               'p_corrected', 'enrich_factor', 'depth', 'study_count', 'pop_count', 'seq_list',
               'neg_log10p_uncorrected', 'neg_log10p_corrected']
    df = df.reindex(columns, axis=1)
    df.to_csv(args.ofile, sep='\t', index=False)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Script for modifying goatools outfile')
    parser.add_argument('-i', action='store', required=True, dest='ifile',
                        help='goatools outfile')
    parser.add_argument('-o', action='store', required=True, dest='ofile',
                        help='result file for api')
    args = parser.parse_args()

    main(args)
