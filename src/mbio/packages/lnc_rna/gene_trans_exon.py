# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, qinjincheng'

from optparse import OptionParser
import os
from collections import defaultdict
import fileinput
import regex

parser = OptionParser(description='Count transcripts number of each gene and exons number of each transcript')
parser.add_option('-i', '--input', dest='input', help='Input directory containing specified GTF files')
parser.add_option('-m', '--method', dest='method', help='Method used for merging GTF files')
parser.add_option('-o', '--output', dest='output', help='Output directory for results')
(opts, args) = parser.parse_args()

def main(dir_in, method, dir_out):
    for gtf, gene_trans, trans_exon in transfer(dir_in, dir_out):
        print 'INFO: start counting transcripts and exons under {} method'.format(method)
        producer(gtf, method, gene_trans, trans_exon)
        if os.path.isfile(gene_trans) and os.path.isfile(trans_exon):
            print 'INFO: succeed in counting transcripts from {}'.format(gene_trans)
            print 'INFO: succeed in counting exons from {}'.format(trans_exon)
        else:
            raise Exception('ERROR: fail to count transcripts and exons from {} or {}'.format(gene_trans, trans_exon))
    else:
        print 'INFO: succeed in counting all transcripts and exons'

def transfer(dir_in, dir_out):
    for file_in in os.listdir(dir_in):
        if file_in in ['old_genes.gtf', 'old_transcripts.gtf', 'new_genes.gtf', 'new_transcripts.gtf']:
            gtf = os.path.join(dir_in, file_in)
            gene_trans = os.path.join(dir_out, '{}.trans'.format(file_in))
            trans_exon = os.path.join(dir_out, '{}.exon'.format(file_in))
            yield gtf, gene_trans, trans_exon

def producer(gtf, method, gene_trans, trans_exon):
    fw1 = open(gene_trans, 'w+')
    fw2 = open(trans_exon, 'w+')
    gene_set_dic = defaultdict(set)
    gene_dic = defaultdict(int)
    trans_dic = defaultdict(int)
    for line in fileinput.input(gtf):
        m = regex.search(r'^[^#]\S*\t(\S+\t){7}.*?gene_id\s+\"(\S+)\";.*\s+transcript_id\s+\"(\S+)\";*', line)
        if m:
            seq_type = m.captures(1)[1]
            txpt_id = m.captures(3)[0]
            gene_id = m.captures(2)[0]
            if str(method).lower() == 'stringtie':
                if seq_type.split("\t")[0] == 'exon':
                    trans_dic[txpt_id] += 1
            elif str(method).lower() == 'cufflinks':
                trans_dic[txpt_id] += 1
            gene_set_dic[gene_id].add(txpt_id)
    for key in gene_set_dic.keys():
        gene_dic[key] = len(gene_set_dic[key])
        fw1.write('{}\t{}\n'.format(key, gene_dic[key]))
    for key in trans_dic.keys():
        fw2.write('{}\t{}\n'.format(key, trans_dic[key]))
    fw1.close()
    fw2.close()

if __name__ == '__main__':
    if opts.input and opts.method and opts.output:
        main(opts.input, opts.method, opts.output)
    else:
        parser.print_help()