# -*- coding : utf-8 -*-
# __author__ : shicaiping
# __date__ : 20200311

import logging
import sqlite3
import os

def main(args):
    if not os.path.exists(args.seq_db):
        logging.error("seq_db file does not exist, please check.")
    if not os.path.exists(args.trans2gene):
        logging.error("transgene file does not exist, please check.")
    t2g = trans2gene(args.trans2gene)
    seq_list = list()
    if args.type == 'gene':
        if args.category == "ref":
            for key in t2g:
                if not t2g[key].startswith(('MSTRG', 'TCONS', 'XLOC')):
                    seq_list.append(t2g[key])
        elif args.category == "new":
            for key in t2g:
                if t2g[key].startswith(('MSTRG', 'TCONS', 'XLOC')):
                    seq_list.append(t2g[key])
        else:
            seq_list = list(set(t2g.values()))
    else:
        if args.category == "ref":
            for key in t2g:
                if not key.startswith(('MSTRG', 'TCONS', 'XLOC')):
                    seq_list.append(key)
        elif args.category == "new":
            for key in t2g:
                if key.startswith(('MSTRG', 'TCONS', 'XLOC')):
                    seq_list.append(key)
        else:
            seq_list = t2g.keys()
    query_seq = parse_db(args.type, seq_list, args.seq_db)
    if 'output_dir' in args and args.output_dir:
        dir = args.output_dir
    else:
        dir = os.getcwd()
    if 'output_name' in args and args.output_name:
        name = args.output_name
    else:
        name = os.path.basename(args.seq_db) + "." + args.type + ".fa"
    output = os.path.join(dir, name)
    with open(output, "w") as w:
        for key in sorted(query_seq):
            w.write(">" + key + "\n")
            start=0
            end=args.width
            seq = query_seq[key]
            while end < len(seq):
                w.write(str(seq[start:end])+"\n")
                start += args.width
                end += args.width
            while end >= len(seq):
                w.write(str(seq[start:])+"\n")
                break

def trans2gene(file):
    t2g = dict()
    with open(file, "r") as f:
        for line in f:
            trans_id = line.strip().split()[0]
            gene_id = line.strip().split()[1]
            if trans_id not in t2g:
                t2g[trans_id] = gene_id
    return t2g

def parse_db(seq_type, seq_id_list, seq_db):
    conn = sqlite3.connect(seq_db)
    cursor = conn.cursor()
    if seq_type == "gene":
        gene_seq = dict()
        query_table = 'gene_seq'
        for seq_id in seq_id_list:
            cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'seq_id', seq_id))
            seq = cursor.fetchall()
            if seq:
                sequence = seq[0][0]
            else:
                sequence = 'None'
            gene_seq[seq_id] = sequence
        return gene_seq
    elif seq_type == "transcript":
        trans_seq = dict()
        query_table = 'trans_seq'
        for seq_id in seq_id_list:
            cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'seq_id', seq_id))
            seq = cursor.fetchall()
            if seq:
                sequence = seq[0][0]
            else:
                sequence = 'None'
            trans_seq[seq_id] = sequence
        return trans_seq
    elif seq_type == "cds":
        cds_seq = dict()
        for seq_id in seq_id_list:
            cursor.execute("SELECT * FROM {} WHERE {}='{}'".format('trans_annot', 'transcript_id', seq_id))
            cds_info = cursor.fetchall()
            fileds = ['transcript_id', 'cds_id', 'pep_id', 'cds_seq', 'pep_seq', 'orf_type']
            if cds_info:
                i = 0
                for each in cds_info:
                    i += 1
                    tmp = dict(zip(fileds, each))
                    cds_id = seq_id + "_cds" + str(i) + " " + tmp['orf_type']
                    cds_sequence = tmp['cds_seq']
                    cds_seq[cds_id] = cds_sequence
        return cds_seq
    elif seq_type == "pep":
        pep_seq = dict()
        for seq_id in seq_id_list:
            cursor.execute("SELECT * FROM {} WHERE {}='{}'".format('trans_annot', 'transcript_id', seq_id))
            pep_info = cursor.fetchall()
            fileds = ['transcript_id', 'cds_id', 'pep_id', 'cds_seq', 'pep_seq', 'orf_type']
            if pep_info:
                i = 0
                for each in pep_info:
                    i += 1
                    tmp = dict(zip(fileds, each))
                    pep_id = seq_id + "_pep" + str(i) + " " + tmp['orf_type']
                    pep_sequence = tmp['pep_seq']
                    pep_seq[pep_id] = pep_sequence
        return pep_seq
    cursor.close()


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description = "This script is to get seqs from seq_db.")
    parser.add_argument('-d', '--seq_db', type=str, required=True, help="sqlite3 build database file")
    parser.add_argument('-t', '--type', type=str, required=True, help="type of sequence to get", choices=['gene', 'transcript', 'cds', 'pep'])
    parser.add_argument('-c', '--category', type=str, required=True, help="type of sequence to get", choices=['ref', 'new', 'all'])
    parser.add_argument('-i', '--trans2gene', type=str, required=True, help="transcript id and gene id correspondence file")
    parser.add_argument('-o', '--output_dir', type=str, help='output directory, default is current dir.')
    parser.add_argument('-n', '--output_name', type=str, help='output file name')
    parser.add_argument('-w', '--width', type=int, default=60, help='width of fasta file per line for output')

    args = parser.parse_args()

    main(args)