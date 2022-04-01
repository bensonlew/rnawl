# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

import os
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd


def parse_blast(table, aln_len, identity):
    # 读取打环序列与线性骨架序列blast表，根据aln_len，identity初步筛选
    # index=["query", "ref", "identity", "aln_len", "mis", "gap", "qstart", "qend", "rstart", "rend", "evalue", "score"]
    data = pd.read_table(table, header=None, names=["query", "ref", "identity", "aln_len", "mis", "gap", "qstart", "qend", "rstart", "rend", "evalue", "score"])
    blast_data = data.loc[(data["aln_len"]>=aln_len) & (data["identity"]>=identity)]
    return blast_data

def parse_log(log_file):
    #根据log文件，如果是环状写入log文件，线状获取比较信息
    genome_info = {}
    my_output_log_record = pd.DataFrame()
    with open(log_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split("\t")
            if line[3] == "True":
                my_output_log_record = my_output_log_record.append({"loc": line[0], "seq": line[1], "from": line[2], "suggest_circle": line[3]}, ignore_index=True)
                print(my_output_log_record)
            if genome_info.has_key(line[0]):
                genome_info[line[1]].append(line[0])
            else:
                genome_info[line[1]] = [line[0]]
    return genome_info, my_output_log_record
def parse(fasta, query, blast_table, log_file, output_prefix, aln_len=1000, identity=95):
    #根据blast筛选的结果，进行搭环处理
    global my_output_log_record
    my_output_fa_record= []
    my_output_log_record = pd.DataFrame()
    seqid_list = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        seqid_list.append(seq_record.id)
    seq_records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    query_records = SeqIO.to_dict(SeqIO.parse(query, "fasta"))
    blast_data = parse_blast(blast_table, aln_len, identity)
    blast_data1 = blast_data.loc[blast_data['rstart'] <= blast_data['rend']]
    blast_data2	= blast_data.loc[blast_data['rstart'] > blast_data['rend']]
    blast_data2.columns =["query", "ref", "identity", "aln_len", "mis", "gap", "qend", "qstart", "rend", "rstart", "evalue", "score"]
    blast_data= pd.concat([blast_data1, blast_data2], ignore_index = True)
    genome_info, my_output_log_record1= parse_log(log_file)
    print(type(my_output_log_record1))
    my_output_log_record= my_output_log_record.append(my_output_log_record1)
    print("aaabbb")
    print(my_output_log_record)
    for seq_id in seqid_list:
        print(seq_id)
        result = blast_data.loc[blast_data["ref"]==seq_id]
        if len(result['aln_len'])<1:
            my_output_fa_record.append(seq_records[seq_id])
            my_output_log_record = my_output_log_record.append(
                {"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "False"},
                ignore_index=True)
            continue
        r_len = len(seq_records[seq_id])
        r_left = result.loc[result['rstart'] <=10]
        aln_len1_max = r_left["aln_len"].max()
        r_left = r_left.loc[r_left['aln_len'] == aln_len1_max]
        if len(r_left['aln_len']) >1:
            identity= r_left["identity"].max()
            r_left2 = r_left.loc[r_left['identity'] == identity]
        elif len(r_left['aln_len']) ==1:
            r_left2 =r_left
        else:
            my_output_fa_record.append(seq_records[seq_id])
            my_output_log_record = my_output_log_record.append(
                {"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "False"},
                ignore_index=True)
            continue
        r_right = result.loc[result['rend'] >= r_len-9]
        aln_len2_max = r_right["aln_len"].max()
        r_right = r_right.loc[r_right['aln_len'] == aln_len2_max]
        if len(r_right['aln_len']) >1:
            identity= r_right["identity"].max()
            r_right2 = r_right.loc[r_right['identity'] == identity]
        elif len(r_right['aln_len']) ==1:
            r_right2 =r_right
        else:
            my_output_fa_record.append(seq_records[seq_id])
            my_output_log_record = my_output_log_record.append(
                {"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "False"},
                ignore_index=True)
            continue
        l_strand = ''
        r_strand = ''
        print(r_left2[0:1])
        if len(r_left2['aln_len']) >1:
            r_left2 = r_left2[0:1]
        if len(r_right2['aln_len']) > 1:
            r_right2 = r_right2[0:1]
        r_left2.set_axis(['a'], inplace=True)
        r_right2.set_axis(['a'], inplace=True)
        if int(r_left2.ix['a', 'rend']) > int(r_left2.ix['a', 'rstart']) and int(r_left2.ix['a', 'qend']) > int(r_left2.ix['a', 'qstart']):
            l_strand = "+"
        elif int(r_left2.ix['a', 'rend']) > int(r_left2.ix['a', 'rstart']) and int(r_left2.ix['a', 'qend']) < int(r_left2.ix['a', 'qstart']):
            l_strand = "-"
        if int(r_right2.ix['a', 'rend']) > int(r_right2.ix['a', 'rstart']) and int(r_right2.ix['a', 'qend']) > int(r_right2.ix['a', 'qstart']):
            r_strand = "+"
        elif int(r_right2.ix['a', 'rend']) > int(r_right2.ix['a', 'rstart']) and int(r_right2.ix['a', 'qend']) < int(r_right2.ix['a', 'qstart']):
            r_strand = "-"
        r_left2_seq_id = r_left2.ix['a', 'query']
        r_right2_seq_id = r_right2.ix['a', 'query']
        if l_strand == "+" and r_strand == "+" and r_left2_seq_id == r_right2_seq_id and r_left2.ix['a', 'rstart'] != r_right2.ix['a', 'rstart'] and r_left2.ix['a', 'qstart'] >= r_right2.ix['a', 'qstart'] and r_left2.ix['a', 'qend'] >= r_right2.ix['a', 'qend']:
            q_start = int(r_left2.ix['a', 'qstart'])
            q_end = int(r_right2.ix['a', 'qend'])
            if q_start <= q_end:
                num =q_end - q_start
                seq = seq_records[seq_id].seq[0:r_len+1-num]
                seq_records[seq_id].seq = seq
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "True"}, ignore_index=True)
            else:
                seq = query_records[r_left2_seq_id].seq[q_end:q_start]
                seq2 = seq_records[seq_id].seq + seq
                seq_records[seq_id].seq = seq2
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "True"}, ignore_index=True)
        elif l_strand == "-" and r_strand == "-" and r_left2_seq_id == r_right2_seq_id and r_left2.ix['a', 'rstart'] != r_right2.ix['a', 'rstart'] and  r_left2.ix['a', 'qend'] <=r_right2.ix['a', 'qend'] and r_left2.ix['a', 'qstart'] <= r_right2.ix['a', 'qstart']:
            print("gaohao")
            q_start = int(r_left2.ix['a', 'qstart'])
            q_end = int(r_right2.ix['a', 'qend'])
            if q_start >= q_end:
                num = q_start - q_end
                seq = seq_records[seq_id].seq[0:r_len + 1 - num]
                seq_records[seq_id].seq =seq
                my_output_fa_record.append(seq_records[seq_id])
                print({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "True"})
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "True"}, ignore_index=True)
            else:
                seq = query_records[r_left2_seq_id].seq[q_start:q_end]
                seq=seq.reverse_complement()
                seq2 =seq_records[seq_id].seq + seq
                seq_records[seq_id].seq=seq2
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "True"}, ignore_index=True)
        else:
            if r_left2.ix['a', 'rstart'] == r_right2.ix['a', 'rstart'] and r_left2.ix['a', 'rend'] == r_right2.ix['a', 'rend']:
                print("aaaa"+str(r_left2))
                print("bbbb"+str(r_right2))
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "False"}, ignore_index=True)
                continue
            df1 = pd.concat([r_left2, r_right2], ignore_index = True)
            aln_len_max = df1["aln_len"].max()
            result = df1.loc[df1["aln_len"] == aln_len_max]
            if len(result['aln_len']) > 1:
                identity = result["identity"].max()
                result2 = result.loc[result['identity'] == identity]
            else:
                result2 = result
            if len(result2['aln_len']) > 1:
                result = result2[0:1]
            result.set_axis(['a'], inplace=True)
            if int(result.ix['a', 'rend']) > int(result.ix['a', 'rstart']) and int(result.ix['a', 'qend']) > int(result.ix['a', 'qstart']) and int(result.ix['a', 'rend']) >= r_len - 9:
                q_start = int(result.ix['a', 'qend'])
                q_seq_id = result.ix['a', 'query']
                seq = query_records[q_seq_id].seq[q_start:]
                seq2 = seq_records[seq_id].seq + seq
                seq_records[seq_id].seq = seq2
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "False"}, ignore_index=True)
            elif result.ix['a', 'rend'] > result.ix['a', 'rstart'] and result.ix['a', 'qend'] > result.ix['a', 'qstart'] and result.ix['a', 'rstart'] <= 10:
                q_end = result.ix['a', 'qstart']
                q_seq_id = result.ix['a', 'query']
                seq = query_records[q_seq_id].seq[:q_end-1]
                seq2 = seq + seq_records[seq_id].seq
                seq_records[seq_id].seq = seq2
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "False"}, ignore_index=True)
            elif result.ix['a', 'rend'] > result.ix['a', 'rstart'] and result.ix['a', 'qend'] < result.ix['a', 'qstart'] and result.ix['a', 'rstart'] <= 10:
                q_end = result.ix['a', 'qstart']
                q_seq_id = result.ix['a', 'query']
                seq = query_records[q_seq_id].seq[q_end:]
                seq = seq.reverse_complement()
                seq2 = seq + seq_records[seq_id].seq
                seq_records[seq_id].seq = seq2
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "False"}, ignore_index=True)
            elif result.ix['a', 'rend'] > result.ix['a', 'rstart'] and result.ix['a', 'qend'] < result.ix['a', 'qstart'] and result.ix['a', 'rend'] >= r_len - 9:
                print("gaohao" + str(result.ix['a', 'qend']))
                q_end = result.ix['a', 'qend']
                q_seq_id = result.ix['a', 'query']
                seq = query_records[q_seq_id].seq[:q_end-1]
                seq = seq.reverse_complement()
                seq2 = seq_records[seq_id].seq + seq
                seq_records[seq_id].seq = seq2
                my_output_fa_record.append(seq_records[seq_id])
                my_output_log_record = my_output_log_record.append({"loc": genome_info[seq_id][0], "seq": seq_id, "from": seq_id, "suggest_circle": "False"}, ignore_index=True)
    print(my_output_log_record)
    my_output_log_record.reindex(columns=["loc", "seq", "from", "suggest_circle"]).to_csv(output_prefix + ".log", sep="\t", header=False, index=False)
    # 生成fasta
    SeqIO.write(my_output_fa_record, output_prefix + ".fa", "fasta")

def _main():
    parser = argparse.ArgumentParser(description='filter fa by json file')
    parser.add_argument('-fa', '--fa', help="fasta file")
    parser.add_argument('-query', '--query', help="query fasta file")
    parser.add_argument('-blast', '--blast', help="blast table in format 6")
    parser.add_argument('-log', '--log', help="log table")
    parser.add_argument('-aln_len', '--aln_len', help="aln_len", type=int)
    parser.add_argument('-identity', '--identity', help="identity", type=float)
    parser.add_argument('-o_prefix', '--o_prefix', help="output prefix, output a new fasta file and a log file")
    args = parser.parse_args()
    parse(args.fa, args.query, args.blast, args.log, args.o_prefix, args.aln_len, args.identity)

if __name__ == "__main__":
    _main()
    '''
        python gap_fill2.py -fa test.fa -query test2.fa -log pre.log -o_prefix test
    '''
