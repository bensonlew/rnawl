# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import pandas as pd
from Bio import SeqIO
import os

parser = OptionParser(description='Select miRNA precursor candidate from lncRNA')
parser.add_option('-i', '--input', dest='input', help='input all lncRNA FASTA file')
parser.add_option('-t', '--tblout', dest='tblout', help='input blastn outfmt6 file')
parser.add_option('-p', '--hairpin', dest='hairpin', help='input miRNA precursor FASTA file')
parser.add_option('-o', '--output', dest='output', help='output candidate lncRNA FASTA file')
(opts, args) = parser.parse_args()

def main(tblout_in, hairpin_fa, file_in, file_out):
    print 'INFO: start processing {}'.format(tblout_in)
    tblout_df = pd.read_table(tblout_in, header=None)
    tblout_df = tblout_df[tblout_df.iloc[:,2] >= 90].copy()
    tblout_df = tblout_df[tblout_df.iloc[:,5] == 0].copy()
    name_map_dict = dict([
        (0, 'query_id'),
        (1, 'subject_id'),
        (2, 'identity'),
        (3, 'alignment_length'),
        (4, 'mismatches'),
        (5, 'gap_openings'),
        (6, 'q_start'),
        (7, 'q_end'),
        (8, 's_start'),
        (9, 's_end'),
        (10, 'e_value'),
        (11, 'bit_score'),
    ])
    tblout_df.rename(columns=name_map_dict, inplace=True)
    print 'INFO: start processing {}'.format(hairpin_fa)
    mir_len_df = pd.DataFrame(
        data=[[seq_record.id, len(seq_record.seq)] for seq_record in SeqIO.parse(hairpin_fa, 'fasta')],
        columns=['subject_id', 'precursor_length']
    )
    print 'INFO: start selecting and writing'
    df = pd.merge(tblout_df, mir_len_df)
    df['length_percent'] = df['alignment_length'] / df['precursor_length']
    df = df[df.loc[:, 'length_percent'] > 0.9].copy()
    SeqIO.write([
        seq_record for seq_record in SeqIO.parse(file_in, 'fasta') if seq_record.id in df['query_id'].values
    ], file_out, 'fasta')
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

if __name__ == '__main__':
    if opts.tblout and opts.hairpin and opts.input and opts.output:
        main(opts.tblout, opts.hairpin, opts.input, opts.output)
    else:
        parser.print_help()