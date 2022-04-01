# coding=utf-8

from Bio import SeqIO
import sqlite3
import os

def extract_transcript_seqs_M(seq, gene_list, output_dir):
    with open(os.path.join(output_dir, "extract_trans_seqs"), "w") as f:
        seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(seq, 'fasta')}
        for i in gene_list:
            try:
                f.write(">{}".format(i) + "\n")
                f.write(seq_dict[i] + "\n")
            except:
                pass

def extract_gene_seqs_G(db, gene_list, output_dir):
    with open(os.path.join(output_dir, "extract_gene_seqs"), "w") as f:
        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        query_table = 'gene_seq'
        for seq_id in gene_list:
            cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'seq_id', seq_id))
            seq = cursor.fetchall()
            f.write(">{}".format(seq_id) + "\n")
            f.write(str(seq[0][0]) + "\n")
        cursor.close()

def extract_tran_seqs_G(db, gene_list, output_dir, type):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    f = open(os.path.join(output_dir, "extract_trans_seqs"), "w")
    if 'cds' in type:
        c = open(os.path.join(output_dir, "extract_cds_seqs"), "w")
    if 'pep' in type:
        p = open(os.path.join(output_dir, "extract_pep_seqs"), "w")
    query_table = 'trans_seq'
    for seq_id in gene_list:
        cursor.execute("SELECT sequence FROM {} WHERE {}='{}'".format(query_table, 'seq_id', seq_id))
        seq = cursor.fetchall()
        f.write(">{}".format(seq_id) + "\n")
        f.write(str(seq[0][0]) + "\n")
        if 'cds' or 'pep' in type:
            cursor.execute("SELECT * FROM {} WHERE {}='{}'".format('trans_annot', 'transcript_id', seq_id))
            cds_info = cursor.fetchall()
            if cds_info:
                if 'cds' in type:
                    c.write(">{}".format(str(cds_info[0][1])) + "\n")
                    c.write(str(cds_info[0][3]) + "\n")
                if 'pep' in type:
                    if cds_info[0][4] == "None":
                        pass
                    else:
                        p.write(">{}".format(str(cds_info[0][2])) + "\n")
                        p.write(str(cds_info[0][4]) + "\n")
    cursor.close()

def extract_seqs(seq=None, ids=None, level=None, type=None, output_dir=None, t2g=None):
    with open(ids, "r") as g:
        all_gene = [line.strip() for line in g.readlines()]
    if level == 'M':
        extract_transcript_seqs_M(seq, all_gene, output_dir)
    if level == 'G':
        extract_gene_seqs_G(seq, all_gene, output_dir)

        if "transcript" in type:
            g2t_dict = dict()
            all_trans = list()
            with open(t2g, 'r') as t:
                for i in t.readlines():
                    tran_id= i.strip().split('\t')[0]
                    gene_id = i.strip().split('\t')[1]
                    if gene_id not in g2t_dict:
                        g2t_dict[gene_id] = [tran_id]
                    else:
                        g2t_dict[gene_id].append(tran_id)
            for i in all_gene:
                all_trans.extend(g2t_dict[i])
            extract_tran_seqs_G(seq, all_trans, output_dir, type)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-s', type=str, metavar="seq", required=True,
                        help="please input fasta file ")
    parser.add_argument('-i', type=str, metavar="gene_transcript_ids", help="id list for extract ", required=True)
    parser.add_argument('-l', type=str, metavar="level", help="gene or transcript level", required=True)
    parser.add_argument('-o', type=str, metavar="output_dir",default=None, help="default is local dir. Output directory name", required=True)
    parser.add_argument('-t', type=str, metavar="extract_type", default=None, help="extract_types ")
    parser.add_argument('-g', type=str, metavar="g2t", default=None, help="g2t file")
    #
    args = parser.parse_args()
    extract_seqs(seq=args.s, ids=args.i, level=args.l, type=args.t, output_dir=args.o, t2g=args.g)

