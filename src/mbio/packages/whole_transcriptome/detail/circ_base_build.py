# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import json
import logging
import pickle
import re
import sqlite3
import os

import pandas as pd
from Bio import SeqIO

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

get_sequence_dict = lambda fasta: {record.id: str(record.seq) for record in SeqIO.parse(fasta, 'fasta')}


def main(args):
    input_dict = json.load(open(args.json))
    biomart_file = input_dict['biomart_file']
    biomart_type = input_dict['biomart_type']
    circrna_fasta = input_dict['circrna_fasta']

    gene2name, gene2des, trans2des = get_des(biomart_file, biomart_type)

    transcript_detail_pk = input_dict['transcript_detail_pk']
    detail_dict = get_detail_dict(input_dict['circrna_detail'], gene2name, gene2des)
    transcript_detail_data = get_transcript_detail_data(detail_dict)
    pickle.dump(transcript_detail_data, open(transcript_detail_pk, 'w'))
    logging.info('succeed in exporting {}'.format(transcript_detail_pk))

    circ_seq_dict = get_sequence_dict(circrna_fasta)
    seqdownload_circ_seq = os.path.join(os.path.dirname(transcript_detail_pk), "seqdownload_circ_seq")
    with open(seqdownload_circ_seq, "wb") as f:
        pickle.dump(circ_seq_dict, f)


def get_des(des, des_type):
    """
    获取蛋白名或描述信息
    """
    gene2name = dict()
    gene2des = dict()
    trans2des = dict()

    with open(des, "rb") as f:
        # f.readline()
        for line in f.readlines():
            line = line.strip().split('\t')
            gene_id = line[0]
            tran_id = line[1]
            if des_type == "type3":
                symbol = ""
                des = line[3]
            elif des_type == "type2":
                symbol = line[2]
                des = line[5]
            elif des_type == "type1":
                symbol = line[2]
                des = line[7]
            else:
                pass

            gene2name[gene_id] = symbol
            gene2des[gene_id] = des
            trans2des[tran_id] = des
    return gene2name, gene2des, trans2des




def get_detail_dict(known_detail, gene2name, gene2des):
    predict_dict = dict()
    if os.path.exists(known_detail):
        with open(known_detail, 'r') as f:
            header = f.readline()

            for line in f:
                cols = line.strip("\n").split("\t")
                try:
                    host_gene = cols[1]
                    chrom = cols[2]
                    start = cols[4]
                    end = cols[5]
                    strand = cols[3]
                    signal = cols[6]
                    type = cols[7]
                except:
                    chrom = "-"
                    start = "-"
                    end = "-"
                    strand = "-"
                    signal = "-"
                    type = "-"
                try:
                    if cols[8] == "yes":
                        kind = "ref"
                    else:
                        kind = "new"
                except:
                    kind = "new"

                predict_dict[cols[0]] = {'kind': kind, 'transcript_id': cols[0], 'chrom': chrom, 'start': start, 'end': end, 'strand': strand, 'signal': signal, 'type': type, 'host_gene': host_gene, 'gene_name': gene2name.get(host_gene, ""), 'description': gene2des.get(host_gene, "")}

    return predict_dict


def get_transcript_detail_data(detail_dict):
    data = list()
    for transcript_id in detail_dict:
        document = detail_dict[transcript_id]
        data.append(document)
    return data





if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Clean bedtools output FASTA file')
    parser.add_argument('--json', action='store', required=True,
                        help='setting JSON file', metavar='<FILE>', dest='json')
    parser.add_argument('--database', action='store', required=True,
                        help='SQLite database file', metavar='<FILE>', dest='database')

    args = parser.parse_args()

    main(args)
