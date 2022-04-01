#!/usr/bin/env python
#coding:utf-8
#fengyitong
#20181106

import argparse
from collections import defaultdict
import os


parser = argparse.ArgumentParser(
    description="""
small_rna质控完文件的碱基偏好性
""")
parser.add_argument(
    "-known",
    dest="known",
    required=True,
    type=str,
    help="已知smallrna的序列文件")
parser.add_argument(
    "-novel",
    dest="novel",
    required=True,
    type=str,
    help="新smallrna的序列文件")
args = parser.parse_args()

def change_to_per(file):
    xls_r = open(file,'r')
    header = xls_r.readline()
    info = xls_r.readlines()
    xls_r.close()
    with open(file,'w') as xls_w, open(file.split('.xls')[0] + '_per.xls','w') as per:
        per.write(header)
        xls_w.write(header)
        for line in info:
            line = line.strip().split('\t')
            id = line[0]
            before = [int(x) for x in line[1:]]
            if sum(before) != 0:
                after = [str(x/(sum(before)+0.0001)) for x in before]
                after = ['1' if float(x) > 0.999 else x for x in after]
                per.write(id + '\t' + '\t'.join(after) + '\n')
                xls_w.write('\t'.join(line) + '\n')
            else:
                continue
    return

def get_first_and_loc_dict(file):
    first_atcg = defaultdict(dict)
    loc_atcg = defaultdict(dict)
    if os.path.getsize(file) == 0:
        return first_atcg, loc_atcg
    else:
        with open(file, 'r') as fasta:
            mirnas = fasta.read().lstrip('>').split('\n>')
            for mirna in mirnas:
                seq = ''.join(mirna.split('\n')[1:])
                seq = seq.replace("T", "U")
                seq = seq.replace("t", "U")
                if not len(seq) in first_atcg:
                    first_atcg[len(seq)] = dict()
                if not seq[0].upper() in first_atcg[len(seq)]:
                    first_atcg[len(seq)][seq[0].upper()] = 0
                first_atcg[len(seq)][seq[0].upper()] += 1
                for i in range(len(seq)):
                    if i not in loc_atcg:
                        loc_atcg[i] = dict()
                    if seq[i].upper() not in loc_atcg[i]:
                        loc_atcg[i][seq[i].upper()] = 0
                    loc_atcg[i][seq[i].upper()] += 1
        return first_atcg, loc_atcg

def creat_loc_and_first(first_atcg, loc_atcg, loc_file, first_file):
    with open(loc_file, 'w') as loc_w, open(first_file, 'w') as length_w:
        loc_w.write('Location\tA\tG\tC\tU\n')
        length_w.write('Length(nt)\tA\tG\tC\tU\n')
        range_f = sorted(first_atcg.keys())
        for i in range_f:
            line = str(i) + '\t'
            for nt in orders:
                if i in first_atcg and nt in first_atcg[i]:
                    line = line + str(first_atcg[i][nt]) + '\t'
                else:
                    line = line + '0\t'
            line = line.strip('\t') + '\n'
            length_w.write(line)
        range_l = sorted(loc_atcg.keys())
        for i in range_l:
            line = str(i + 1) + '\t'
            for nt in orders:
                if i in loc_atcg and nt in loc_atcg[i]:
                    line = line + str(loc_atcg[i][nt]) + '\t'
                else:
                    line = line + '0\t'
            line = line.strip('\t') + '\n'
            loc_w.write(line)
    change_to_per(loc_file)
    change_to_per(first_file)

def merge_some_dict(know_dict, novel_dict):
    merge_dict = defaultdict(dict)
    merge_keys = set(know_dict.keys() + novel_dict.keys())
    for key in merge_keys:
        merge_dict[key] = dict()
        for nt in orders:
            merge_dict[key][nt] = 0
            know_num = 0
            if key in know_dict and nt in know_dict[key]:
                know_num = know_dict[key][nt]
            novel_num = 0
            if key in novel_dict and nt in novel_dict[key]:
                novel_num = novel_dict[key][nt]
            merge_dict[key][nt] = know_num + novel_num
    return merge_dict


orders=["A","G","C","U"]
first_atcg_known, loc_atcg_known = get_first_and_loc_dict(args.known)
first_atcg_novel, loc_atcg_novel = get_first_and_loc_dict(args.novel)
first_atcg_all = merge_some_dict(first_atcg_known, first_atcg_novel)
loc_atcg_all = merge_some_dict(loc_atcg_known, loc_atcg_novel)
creat_loc_and_first(first_atcg_known, loc_atcg_known, 'known_loc_bias.xls', 'known_first_bias.xls')
creat_loc_and_first(first_atcg_novel, loc_atcg_novel, 'novel_loc_bias.xls', 'novel_first_bias.xls')
creat_loc_and_first(first_atcg_all, loc_atcg_all, 'all_loc_bias.xls', 'all_first_bias.xls')
