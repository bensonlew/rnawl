# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from optparse import OptionParser
import xml.etree.ElementTree as ET
import os
import csv

parser = OptionParser(description='Generate gene BLAST id file')
parser.add_option('--t2g', dest='t2g', help='input T2G file')
parser.add_option('--tids', dest='tids', type=str, help='input transcript ids file')
parser.add_option('--output_t', dest='output_t', type=str, help='output tran ids  file')
parser.add_option('--output_g', dest='output_g', type=str, help='output gene ids  file')
parser.add_option('--output_format', dest='output_format', type=str, help='output gene ids  file')
parser.add_option('--join', dest='join', default='; ', type=str, help='output gene ids  file')

(opts, args) = parser.parse_args()

head = '''<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
'''

def main(t2g_tsv, tids, tids_out, gids_out, aformat, join_str):
    t2g_dict = dict(line.strip().split('\t')[:2] for line in open(t2g_tsv))

    tran_dict = dict()
    gene_dict = dict()
    with open(tids, 'r') as f:
        for id_dict in csv.DictReader(f, delimiter='\t'):
            tran_id = id_dict['transcript_id']
            if "nr" in aformat:
                if id_dict['nr'] != "":
                    id_des = aformat.format(nr=id_dict['nr'], description=id_dict['description'])
            elif 'swissprot' in aformat:
                if id_dict['swissprot'] != "":
                    id_des = aformat.format(swissprot=id_dict['swissprot'], description=id_dict['description'])
            elif 'uniprot' in aformat:
                if id_dict['uniprot'] != "":
                    id_des = aformat.format(uniprot=id_dict['uniprot'], description=id_dict['description'])

            # print tran_id
            if tran_id in tran_dict:
                tran_dict[tran_id].add(id_des)
            else:
                tran_dict[tran_id] = set([id_des])

            if tran_id in t2g_dict:
                gene_id = t2g_dict[tran_id]
                if gene_id in gene_dict:
                    gene_dict[gene_id].add(id_des)
                else:
                    gene_dict[gene_id] = set([id_des])

    # print "read nr ids"
    with open(tids_out, 'w') as f:
        for k, v in tran_dict.items():
            gos = join_str.join(list(v))
            f.write("{}\t{}\n".format(k, gos))

    with open(gids_out, 'w') as f:
        for k, v in gene_dict.items():
            gos = join_str.join(list(v))
            f.write("{}\t{}\n".format(k, gos))


if __name__ == '__main__':
    for i in ['t2g', 'tids', 'output_g', 'output_t', 'output_format']:
        print i, getattr(opts, i)
    if all(map(hasattr, [opts] * 5, ['t2g', 'tids', 'output_g', 'output_t', 'output_format'])):
        main(opts.t2g, opts.tids, opts.output_t, opts.output_g, opts.output_format, opts.join)
    else:
        parser.print_help()
