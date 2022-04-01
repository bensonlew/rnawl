# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
from Bio import SeqIO
import os

parser = OptionParser(description='Process BLAST tblout')
parser.add_option('-i', '--input', dest='input', help='input BLAST tblout')
parser.add_option('-k', '--known', dest='known', help='input known lncRNA id list file')
parser.add_option('-n', '--novel', dest='novel', help='input novel lncRNA id list file')
parser.add_option('-t', '--t2g', dest='t2g', help='input TXPT2GENE file')
parser.add_option('-p', '--hairpin', dest='hairpin', help='input miRNA precursor FASTA file')
parser.add_option('-o', '--output', dest='output', help='output tabular file')
(opts, args) = parser.parse_args()

def main(file_in, known_id_file, novel_id_file, t2g_file, hairpin_fa, file_out):
    print 'INFO: start reading {}'.format(file_in)
    columns = [
        'lncRNA_id', 'pre-miRNA_name', 'Identity', 'Alignment_length', 'mismatch', 'gap_openings',
        'lncRNA_start', 'lncRNA_end', 'pre-miRNA_start', 'pre-miRNA_end', 'E-value', 'Score'
    ]
    df = pd.DataFrame(
        data=[line.strip().split('\t') for line in open(file_in)],
        columns=columns
    )
    df['Identity'] = df['Identity'].apply(float)
    df = df[df['Identity'] >= 90]
    df['gap_openings'] = df['gap_openings'].apply(int)
    df = df[df['gap_openings'] == 0]
    df['E-value'] = df['E-value'].apply(float)
    df = df.sort_values(by=['E-value'])
    print 'INFO: start reading {} and {}'.format(known_id_file, novel_id_file)
    df = df.merge(pd.concat([
        pd.DataFrame(
            data=[(i, 'known') for i in [line.strip() for line in open(known_id_file)]],
            columns=['lncRNA_id', 'lncRNA_type']
        ),
        pd.DataFrame(
            data=[(i, 'novel') for i in [line.strip() for line in open(novel_id_file)]],
            columns=['lncRNA_id', 'lncRNA_type']
        )
    ], ignore_index=True), on='lncRNA_id', how='left')
    print 'INFO: start reading {}'.format(t2g_file)
    df = df.merge(pd.DataFrame(
        data=[line.strip().split('\t') for line in open(t2g_file)], columns=['lncRNA_id', 'lncRNA_gene_id']
    ), on='lncRNA_id', how='left')
    print 'INFO: start reading {}'.format(hairpin_fa)
    df = df.merge(pd.DataFrame(
        data=[[seq_record.id, len(seq_record)] for seq_record in SeqIO.parse(hairpin_fa, 'fasta')],
        columns=['pre-miRNA_name', 'pre-miRNA_length']
    ), on='pre-miRNA_name', how='left')
    df['Alignment_length'] = df['Alignment_length'].apply(int)
    df['pre-miRNA_length'] = df['pre-miRNA_length'].apply(int)
    df = df[df['Alignment_length'] / df['pre-miRNA_length'] > 0.9]
    df = df.reindex(columns=[
        'lncRNA_id', 'lncRNA_gene_id', 'pre-miRNA_name', 'Identity', 'pre-miRNA_length', 'Alignment_length',
        'E-value', 'Score', 'lncRNA_start', 'lncRNA_end',
        'pre-miRNA_start', 'pre-miRNA_end', 'mismatch', 'lncRNA_type'
    ])
    df = df.rename(columns={
        'lncRNA_id': 'lncrna_id',
        'lncRNA_gene_id': 'lncrna_gene_id',
        'pre-miRNA_name': 'pre_mirna_name',
        'Identity': 'identity',
        'pre-miRNA_length': 'pre_mirna_length',
        'Alignment_length': 'alignment_length',
        'E-value': 'e_value',
        'Score': 'score',
        'lncRNA_start': 'lncrna_start',
        'lncRNA_end': 'lncrna_end',
        'pre-miRNA_start': 'pre_mirna_start',
        'pre-miRNA_end': 'pre_mirna_end',
        'mismatch': 'mismatch',
        'lncRNA_type': 'lncrna_type',
    })
    df.to_csv(file_out, sep='\t', index=False)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.input and opts.known and opts.novel and opts.t2g and opts.hairpin and opts.output:
        main(opts.input, opts.known, opts.novel, opts.t2g, opts.hairpin, opts.output)
    else:
        parser.print_help()