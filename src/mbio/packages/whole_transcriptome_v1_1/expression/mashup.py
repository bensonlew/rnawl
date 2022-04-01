# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def main(args):
    program = args.program
    list_file = args.input
    level = args.level
    output_dir = args.output
    merge_func = {'rsem': merge_rsem_results, 'kallisto': merge_kallisto_results, 'salmon': merge_salmon_results}[
        program]
    merge_func(list_file, level, output_dir)


def merge_rsem_results(list_file, level, output_dir):
    tpm_dfs = list()
    fpkm_dfs = list()
    count_dfs = list()
    tpm_matrix_file = os.path.join(output_dir, '{}.tpm.txt'.format(level))
    fpkm_matrix_file = os.path.join(output_dir, '{}.fpkm.txt'.format(level))
    count_matrix_file = os.path.join(output_dir, '{}.count.txt'.format(level))
    seq_id = 'transcript_id' if level == 'T' else 'gene_id'
    for line in open(list_file):
        sample, filepath = line.strip().split('\t')
        quant_df = pd.read_table(filepath)
        tpm_df = quant_df.reindex([seq_id, 'TPM'], axis=1)
        tpm_dfs.append(tpm_df.rename({seq_id: 'seq_id', 'TPM': sample}, axis=1).set_index('seq_id'))
        fpkm_df = quant_df.reindex([seq_id, 'FPKM'], axis=1)
        fpkm_dfs.append(fpkm_df.rename({seq_id: 'seq_id', 'FPKM': sample}, axis=1).set_index('seq_id'))
        count_df = quant_df.reindex([seq_id, 'expected_count'], axis=1)
        count_dfs.append(count_df.rename({seq_id: 'seq_id', 'expected_count': sample}, axis=1).set_index('seq_id'))
    else:
        tpm_matrix_df = pd.concat(tpm_dfs, axis=1)
        tpm_matrix_df.index.name = 'seq_id'
        tpm_matrix_df.to_csv(tpm_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(tpm_matrix_file))
        fpkm_matrix_df = pd.concat(fpkm_dfs, axis=1)
        fpkm_matrix_df.index.name = 'seq_id'
        fpkm_matrix_df.to_csv(fpkm_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(fpkm_matrix_file))
        count_matrix_df = pd.concat(count_dfs, axis=1)
        count_matrix_df.index.name = 'seq_id'
        count_matrix_df.to_csv(count_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(count_matrix_file))


def merge_kallisto_results(list_file, level, output_dir):
    tpm_dfs = list()
    fpkm_dfs = list()
    count_dfs = list()
    tpm_matrix_file = os.path.join(output_dir, '{}.tpm.txt'.format(level))
    fpkm_matrix_file = os.path.join(output_dir, '{}.fpkm.txt'.format(level))
    count_matrix_file = os.path.join(output_dir, '{}.count.txt'.format(level))
    seq_id = 'transcript_id' if level == 'T' else 'gene_id'
    for line in open(list_file):
        sample, filepath = line.strip().split('\t')
        quant_df = pd.read_table(filepath)
        tpm_df = quant_df.reindex([seq_id, 'tpm'], axis=1)
        tpm_dfs.append(tpm_df.rename({seq_id: 'seq_id', 'tpm': sample}, axis=1).set_index('seq_id'))
        fpkm_df = quant_df.reindex([seq_id, 'fpkm'], axis=1)
        fpkm_dfs.append(fpkm_df.rename({seq_id: 'seq_id', 'fpkm': sample}, axis=1).set_index('seq_id'))
        count_df = quant_df.reindex([seq_id, 'count'], axis=1)
        count_dfs.append(count_df.rename({seq_id: 'seq_id', 'count': sample}, axis=1).set_index('seq_id'))
    else:
        tpm_matrix_df = pd.concat(tpm_dfs, axis=1)
        tpm_matrix_df.index.name = 'seq_id'
        tpm_matrix_df.to_csv(tpm_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(tpm_matrix_file))
        fpkm_matrix_df = pd.concat(fpkm_dfs, axis=1)
        fpkm_matrix_df.index.name = 'seq_id'
        fpkm_matrix_df.to_csv(fpkm_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(fpkm_matrix_file))
        count_matrix_df = pd.concat(count_dfs, axis=1)
        count_matrix_df.index.name = 'seq_id'
        count_matrix_df.to_csv(count_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(count_matrix_file))


def merge_salmon_results(list_file, level, output_dir):
    tpm_dfs = list()
    fpkm_dfs = list()
    count_dfs = list()
    tpm_matrix_file = os.path.join(output_dir, '{}.tpm.txt'.format(level))
    fpkm_matrix_file = os.path.join(output_dir, '{}.fpkm.txt'.format(level))
    count_matrix_file = os.path.join(output_dir, '{}.count.txt'.format(level))
    for line in open(list_file):
        sample, filepath = line.strip().split('\t')
        quant_df = pd.read_table(filepath)
        # method for computing TPM
        # denom = sum(quant_df['NumReads'] / quant_df['EffectiveLength'])
        # quant_df['TPM'] = quant_df['NumReads'] / quant_df['EffectiveLength'] / denom * 1e6
        # method for computing FPKM
        denom = quant_df['EffectiveLength'].sum()
        quant_df['FPKM'] = quant_df['NumReads'] / quant_df['EffectiveLength'] / denom * 1e9
        tpm_df = quant_df.reindex(['Name', 'TPM'], axis=1)
        tpm_dfs.append(tpm_df.rename({'Name': 'seq_id', 'TPM': sample}, axis=1).set_index('seq_id'))
        fpkm_df = quant_df.reindex(['Name', 'FPKM'], axis=1)
        fpkm_dfs.append(fpkm_df.rename({'Name': 'seq_id', 'FPKM': sample}, axis=1).set_index('seq_id'))
        count_df = quant_df.reindex(['Name', 'NumReads'], axis=1)
        count_dfs.append(count_df.rename({'Name': 'seq_id', 'NumReads': sample}, axis=1).set_index('seq_id'))
    else:
        tpm_matrix_df = pd.concat(tpm_dfs, axis=1)
        tpm_matrix_df.index.name = 'seq_id'
        tpm_matrix_df.to_csv(tpm_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(tpm_matrix_file))
        fpkm_matrix_df = pd.concat(fpkm_dfs, axis=1)
        fpkm_matrix_df.index.name = 'seq_id'
        fpkm_matrix_df.to_csv(fpkm_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(fpkm_matrix_file))
        count_matrix_df = pd.concat(count_dfs, axis=1)
        count_matrix_df.index.name = 'seq_id'
        count_matrix_df.to_csv(count_matrix_file, sep='\t')
        logger.info('succeed in exporting {}'.format(count_matrix_file))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate expression quantity matrix from multiple quantity results')
    parser.add_argument('-p', action='store', choices=['rsem', 'kallisto', 'salmon'], required=True,
                        help='program for quantification', dest='program')
    parser.add_argument('-i', action='store', required=True,
                        help='quantity result LIST file', metavar='<FILE>', dest='input')
    parser.add_argument('-l', action='store', choices=['G', 'T'], required=True,
                        help='expression level', dest='level')
    parser.add_argument('-o', action='store', required=True,
                        help='output directory', metavar='<DIR>', dest='output')

    args = parser.parse_args()

    main(args)
