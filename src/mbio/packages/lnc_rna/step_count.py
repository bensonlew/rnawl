# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, qinjincheng'

from optparse import OptionParser
from Bio import SeqIO
import os

parser = OptionParser(description='Export statistical results of transcripts count')
parser.add_option('-f', '--fasta_in', dest='fasta_in', help='Input FASTA file of merged transcripts')
parser.add_option('-d', '--dir_in', dest='dir_in', help='Input directory for holding seq length info file')
parser.add_option('-g', '--group_num', dest='group_num', help='Group Number of length distribution')
parser.add_option('-s', '--steps', dest='steps', help='Step size of the distribution')
parser.add_option('-o', '--dir_out', dest='dir_out', help='Output directory for statistical results')
(opts, args) = parser.parse_args()

def main(fasta_in, dir_in, group_num, steps, dir_out):
    seq_len_file = os.path.join(dir_in, '{}.txt'.format(os.path.basename(fasta_in)))
    ret = seq_len_stat(fasta_in, seq_len_file)
    if ret:
        for step in steps.split(','):
            stat_out = os.path.join(dir_out, 'trans_count_stat_{}.txt'.format(step))
            ret = step_count(fasta_in, seq_len_file, int(group_num), int(step), stat_out)
            if not ret:
                raise Exception('ERROR: fail to create {}'.format(stat_out))
        else:
            print 'INFO: finish exporting statistical result of transcripts counts'
    else:
        raise Exception('ERROR: fail to create {}'.format(seq_len_file))

def seq_len_stat(fasta_in, seq_len_out):
    print 'INFO: start recording seq length from {}'.format(fasta_in)
    with open(seq_len_out, 'w') as f:
        for seq_record in SeqIO.parse(fasta_in, 'fasta'):
            seq_id = seq_record.description.strip().split(' ')[0]
            f.write('{}\t{}\n'.format(seq_id, len(seq_record)))
    if os.path.getsize(seq_len_out) > 0:
        print 'INFO: succeed in creating {}'.format(seq_len_out)
        return True

def step_count(fasta_in, seq_len_in, group_num, step, stat_out):
    print 'INFO: start counting the length distribution by {} step at {} groups'.format(step, group_num)
    with open(seq_len_in) as r, open(stat_out, 'a') as w:
        sample_name = os.path.basename(fasta_in).split('.fa')[0]
        w.write('{}\n'.format(sample_name))
        trans_list = list()
        amount_group = list()
        element_set = set()
        for line in r:
            line = line.strip().split('\t')
            number = line[1]
            trans_list.append(number)
        for f in trans_list:
            for i in range(group_num):
                if (int(f) >= (i * step)) and (int(f) < ((i+1) * step)):
                    amount_group.append(i)
                element_set.add(i)
        amount_group.sort()
        top_sum = 0
        for i in element_set:
            num_statistics = amount_group.count(i)
            if i == 0:
                w.write('{}~{}\t{}\n'.format(i * step, (i + 1) * step, num_statistics))
                top_sum += num_statistics
            elif i < group_num - 1:
                w.write('{}~{}\t{}\n'.format(i * step + 1, (i + 1) * step, num_statistics))
                top_sum += num_statistics
            else:
                w.write('>{}\t{}\n'.format(i * step, len(trans_list) - top_sum))
                w.write('total\t{}\n'.format(len(trans_list)))
    if os.path.getsize(stat_out) > 0:
        print 'INFO: succeed in creating {}'.format(stat_out)
        return True

if __name__ == '__main__':
    if opts.fasta_in and opts.dir_in and opts.group_num and opts.steps and opts.dir_out:
        main(opts.fasta_in, opts.dir_in, opts.group_num, opts.steps, opts.dir_out)
    else:
        parser.print_help()