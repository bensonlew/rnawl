# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Conduct statistics of the quantitative results of expression module')
parser.add_option('-i', dest='input', type=str, help='input target document location-name relationship file')
parser.add_option('-m', dest='method', choices=['rsem', 'salmon', 'kallisto'], help='quantitative method')
parser.add_option('-t', dest='type', choices=['T', 'G'], help='expression level')
parser.add_option('-c', dest='map', type=str, help='output clean T2G relationship file')
parser.add_option('-o', dest='output', type=str, help='output directory containing expression matrix table')
(opts, args) = parser.parse_args()

def main(file_in, method, exp_type, t2g, dir_out):
    print 'INFO: start reading {}'.format(file_in)
    n2l_dict = dict([line.strip().split('\t')[::-1] for line in open(file_in)])
    if method == 'rsem':
        outputs = build_rsem_matrix(n2l_dict, exp_type, dir_out)
    elif method == 'salmon':
        outputs = build_salmon_matrix(n2l_dict, exp_type, dir_out)
    elif method == 'kallisto':
        outputs = build_kallisto_matrix(n2l_dict, exp_type, t2g, dir_out)
    if outputs:
        for file_out in outputs:
            if os.path.getsize(file_out) > 0:
                print 'INFO: succeed in exporting {}'.format(file_out)
    else:
        raise Exception('ERROR: fail to generate result in {}'.format(dir_out))

def build_rsem_matrix(n2l_dict, exp_type, dir_out):
    dfs = [(name, pd.read_table(n2l_dict[name], usecols=[0, 4, 5, 6], index_col=0)) for name in sorted(n2l_dict.keys())]
    outputs = list()
    for k, v in {'expected_count': 'count', 'TPM': 'tpm', 'FPKM': 'fpkm'}.items():
        output = os.path.join(dir_out, 'rsem.{}.tsv'.format(v))
        matrix = pd.concat([df[k].rename(name) for name, df in dfs], axis=1)
        matrix.index.rename('seq_id', inplace=True)
        matrix.to_csv(output, sep='\t')
        outputs.append(output)
    else:
        return outputs

def build_salmon_matrix(n2l_dict, exp_type, dir_out):
    dfs = [(name, pd.read_table(n2l_dict[name], usecols=[0, 3, 4], index_col=0)) for name in sorted(n2l_dict.keys())]
    outputs = list()
    for k, v in {'TPM': 'tpm', 'NumReads': 'count'}.items():
        output = os.path.join(dir_out, 'salmon.{}.{}.tsv'.format(exp_type, v))
        matrix = pd.concat([df[k].rename(name) for name, df in dfs], axis=1)
        matrix.index.rename('seq_id', inplace=True)
        matrix.to_csv(output, sep='\t')
        outputs.append(output)
    else:
        return outputs

def build_kallisto_matrix(n2l_dict, exp_type, t2g, dir_out):
    dfs = [(name, pd.read_table(n2l_dict[name], usecols=[0, 3, 4], index_col=0)) for name in sorted(n2l_dict.keys())]
    outputs = list()
    if exp_type == 'G':
        mdf = pd.read_table(t2g, index_col=0, names=['target_id', 'gene_id'])
    for k, v in {'est_counts': 'count', 'tpm': 'tpm'}.items():
        output = os.path.join(dir_out, 'kallisto.{}.{}.tsv'.format(exp_type, v))
        if exp_type == 'T':
            matrix = pd.concat([df[k].rename(name) for name, df in dfs], axis=1)
        elif exp_type == 'G':
            matrix = pd.concat([pd.concat([df[k].rename(name) for name, df in dfs], axis=1), mdf], axis=1)
            matrix = matrix.groupby(matrix['gene_id']).sum()
        matrix.index.rename('seq_id', inplace=True)
        matrix.to_csv(output, sep='\t')
        outputs.append(output)
    else:
        return outputs

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 5, ['input', 'method', 'type', 'map', 'output'])):
        main(opts.input, opts.method, opts.type, opts.map, opts.output)
    else:
        parser.print_help()
