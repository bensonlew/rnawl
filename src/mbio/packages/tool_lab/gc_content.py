#!/usr/bin/env python
# coding: utf8

import re
import argparse

# 读取文件，统计长度
def gc_content(fasta_file,out_name):
    dict_GC = {};dict_N = {};dict_AT = {};dict_length = {};all_key = []
    seq_num = 0;all_GC = 0;all_N = 0;all_AT = 0;seqs_sum=0
    with open(fasta_file,'r') as read_fasta, open(out_name +".all.result.xls",'w') as f, open(out_name + ".detail.result.xls",'w') as t:
        f.write("sample name" + "\t" + "seq count" + "\t" + "GC content (%)" + "\t" + "N ratio(%)" + "\t" + "AT/GC ratio" + "\n")
        t.write("sample name" + "\t" + "scaffold name" + "\t" + "length (bp)" + "\t" + "GC content (%)" + "\t" + "N ratio (%)" + "\t" + "AT/GC ratio" + "\n")
        for line in read_fasta:
            if line[0] == '>':
                key = line.split(" ")[0].strip('[ >\n]')
                all_key.append(key)
                dict_GC[key] = 0
                dict_N[key] = 0
                dict_AT[key] = 0
                dict_length[key] = 0
                seqs_sum += 1
            else:
                value = line.strip()
                seqs_len = len(value)
                dict_length[key] += seqs_len
                dict_GC[key] += (value.count('g') + value.count('G') + value.count('c') + value.count('C'))
                dict_N[key] += (value.count('n') + value.count('N'))
                dict_AT[key] += (value.count('a') + value.count('A') + value.count('t') + value.count('T'))
        for i in all_key:
            all_GC += dict_GC[i]
            all_N += dict_N[i]
            all_AT += dict_AT[i]
            t.write(out_name + "\t" + i + "\t" + str(dict_length[i]) + "\t" + str(
                round((float(dict_GC[i]) / dict_length[i]) * 100, 2)) + "\t" + str(round((float(dict_N[i]) / dict_length[i]) * 100, 2)) + "\t" + str(
                round(float(dict_AT[i]) / dict_GC[i], 2)) + "\n")
        all_base = all_GC + all_N + all_AT
        f.write(out_name + "\t" + str(seqs_sum) + "\t" + str(round((float(all_GC) / all_base) * 100, 2)) + "\t" + str(
            round((float(all_N) / all_base) * 100, 2)) + "\t" + str(round((float(all_AT) / all_GC) , 2)) + "\n")

def _main():
    parser = argparse.ArgumentParser(description='GC Content')
    parser.add_argument('-i', '--fasta_file', help="fasta_file")
    parser.add_argument('-o', '--out_name', help="out_name")
    args = parser.parse_args()
    gc_content(args.fasta_file,args.out_name)

if __name__ == "__main__":
    _main()
