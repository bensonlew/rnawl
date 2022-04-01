# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
from optparse import OptionParser

import numpy as np
import pandas as pd

parser = OptionParser(description='Prepare data and setting file for STEM')
parser.add_option('--matrix', dest='matrix', type=str, help='raw data matrix')
parser.add_option('--geneset', dest='geneset', type=str, help='query geneset')
parser.add_option('--pheno', dest='pheno', type=str, help='pheno table')
parser.add_option('--log', dest='log', type=int, help='log base of normalization')
parser.add_option('--method', dest='method', choices=['SCM', 'K'], help='clustering method')
parser.add_option('--number', dest='number', type=int, help='maximum number of model profiles')
parser.add_option('--unit', dest='unit', type=int, help='maximum unit change in model profiles between time points')
parser.add_option('--significance', dest='significance', type=float, help='significance level')
parser.add_option('--correction', dest='correction', choices=['Bonferroni', 'FDR', 'None'], help='correction method')
parser.add_option('--clusters', dest='clusters', type=int, help='number of clusters K')
parser.add_option('--starts', dest='starts', type=int, help='number of random starts')
parser.add_option('--output', dest='output', type=str, help='output setting file')
(opts, args) = parser.parse_args()

MAIN_SETTING = '''#Main Input:
Data_File	{data}
Gene_Annotation_Source	No annotations
Gene_Annotation_File	
Cross_Reference_Source	No cross references
Cross_Reference_File	
Gene_Location_Source	No Gene Locations
Gene_Location_File	
Clustering_Method[STEM Clustering Method,K-means]	{method}
'''

SCM_SETTING = '''Maximum_Number_of_Model_Profiles	{number}
Maximum_Unit_Change_in_Model_Profiles_between_Time_Points	{unit}
Normalize_Data[Log normalize data,Normalize data,No normalization/add 0]	{scale}
Spot_IDs_included_in_the_data_file	false

#Filtering:
Maximum_Number_of_Missing_Values	0
Minimum_Correlation_between_Repeats	0.0
Minimum_Absolute_Log_Ratio_Expression	1
Change_should_be_based_on[Maximum-Minimum,Difference From 0]	Maximum-Minimum
Pre-filtered_Gene_File	

#Model Profiles
Maximum_Correlation	1.0
Number_of_Permutations_per_Gene	50
Maximum_Number_of_Candidate_Model_Profiles	1000000
Significance_Level	{significance}
Correction_Method[Bonferroni,False Discovery Rate,None]	{correction}
Permutation_Test_Should_Permute_Time_Point_0	true

#Clustering Profiles:
Clustering_Minimum_Correlation	0.7
Clustering_Minimum_Correlation_Percentile	0.0
'''

K_SETTING = '''Number_of_Clusters_K	{clusters}
Number_of_Random_Starts	{starts}
Normalize_Data[Log normalize data,Normalize data,No normalization/add 0]	{scale}
Spot_IDs_included_in_the_data_file	false

#Filtering:
Maximum_Number_of_Missing_Values	0
Minimum_Correlation_between_Repeats	0.0
Minimum_Absolute_Log_Ratio_Expression	1
Change_should_be_based_on[Maximum-Minimum,Difference From 0]	Maximum-Minimum
Pre-filtered_Gene_File	
'''

SUB_SETTING = {'SCM': SCM_SETTING, 'K': K_SETTING}


def main(opts):
    out_matrix = os.path.join(os.path.dirname(opts.output), 'data.tsv')
    print 'INFO: start generating {}'.format(out_matrix)
    print 'INFO: opts -> ({})'.format(opts)
    log_base = opts.log if hasattr(opts, 'log') else None
    if hasattr(opts, 'geneset'):
        export_group_matrix(opts.matrix, opts.pheno, log_base, out_matrix, opts.geneset)
    else:
        export_group_matrix(opts.matrix, opts.pheno, log_base, out_matrix)
    attrs = ['method', 'number', 'unit', 'significance', 'correction', 'clusters', 'starts']
    kwargs = dict((attr, getattr(opts, attr)) for attr in attrs if hasattr(opts, attr))
    kwargs.update({'scale': 'Normalize data' if log_base else 'Log normalize data'})
    print 'INFO: start initializing paramters'
    print 'INFO: kwargs -> ({})'.format(kwargs)
    params = initialize_paramters(**kwargs)
    params.update({'data': out_matrix})
    print 'INFO: params -> ({})'.format(params)
    template = MAIN_SETTING + SUB_SETTING[opts.method]
    open(opts.output, 'w').write(template.format(**params))
    if os.path.getsize(opts.output):
        print 'INFO: succeed in exporting {}'.format(opts.output)


def export_group_matrix(raw_matrix, pheno_table, log_base, out_matrix, geneset_list=None):
    df = pd.read_table(raw_matrix)
    phedf = pd.read_table(pheno_table)
    phedf = phedf.sort_values(['Order', 'Group', 'Sample'])
    samples = phedf['Sample'].tolist()
    groups = phedf['Group'].tolist()
    ## 增加判断，兼容全转录组表结构
    if 'seq_id' in df.columns.values:
        df = df.reindex(['seq_id'] + samples, axis=1).set_index('seq_id')
    # elif 'gene_id' in df.columns.values and len(df['gene_id']) == len(set(df['gene_id'])):
    #     df.rename(columns={'gene_id': 'seq_id'}, inplace=True)
    #     df = df.reindex(['seq_id'] + samples, axis=1).set_index('seq_id')
    # elif 'transcript_id' in df.columns.values:
    #     df.rename(columns={'transcript_id': 'seq_id'}, inplace=True)
    #     df = df.reindex(['seq_id'] + samples, axis=1).set_index('seq_id')
    else:
        raise Exception("Can't find seq_id!")
    if geneset_list:
        genesets = [line.strip() for line in open(geneset_list)]
        df = df.query('index in @genesets')
    columns = pd.MultiIndex.from_arrays([groups, samples], names=['group', 'sample'])
    df = pd.DataFrame(np.array(df), index=df.index, columns=columns)
    df = df.groupby(level='group', axis=1).mean()
    group = list()
    for g in groups:
        if g not in group:
            group.append(g)
    df = df.reindex(group, axis=1)
    if log_base:
        df = df.apply(lambda x: np.log(x + 1) / np.log(log_base))
    df.to_csv(out_matrix, sep='\t')


def initialize_paramters(**kwargs):
    if kwargs['method'] == 'SCM':
        corr_map = {'Bonferroni': 'Bonferroni', 'FDR': 'False Discovery Rate', 'None': 'None'}
        parameters = {
            'method': 'STEM Clustering Method',
            'number': kwargs['number'],
            'unit': kwargs['unit'],
            'significance': kwargs['significance'],
            'correction': corr_map[kwargs['correction']],
            'scale': kwargs['scale']
        }
    elif kwargs['method'] == 'K':
        parameters = {
            'method': 'K-means',
            'clusters': kwargs['clusters'],
            'starts': kwargs['starts'],
            'scale': kwargs['scale']
        }
    return parameters


if __name__ == '__main__':
    if all(map(hasattr, [opts] * 4, ['matrix', 'pheno', 'method', 'output'])):
        main(opts)
    else:
        parser.print_help()
