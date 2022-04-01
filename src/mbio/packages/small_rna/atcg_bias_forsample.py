#!/usr/bin/env python
#coding:utf-8
#fengyitong
#20180621

import argparse
import ConfigParser
from collections import defaultdict


parser = argparse.ArgumentParser(
    description="""
small_rna质控完文件的碱基偏好性
""")
parser.add_argument(
    "-fa",
    dest="fasta",
    required=True,
    type=str,
    help="质控完的序列文件")
parser.add_argument(
    "-cfg",
    dest="config",
    required=True,
    type=str,
    help="有命名关系的config文件")
parser.add_argument(
    "-min",
    dest="min",
    type=str,
    default = '18',
    help='最短序列的长度')
parser.add_argument(
    "-max",
    dest="max",
    type=str,
    default = '32',
    help='最长序列的长度')
args = parser.parse_args()

cfg = ConfigParser.ConfigParser()
cfg.read(args.config)
names = cfg.options('NAME')
names = [x.lower() for x in names]

def change_to_per(file):
    with open(file,'r') as xls_r, open(file.strip('.xls') + '_per.xls','w') as per:
        header = xls_r.readline()
        per.write(header)
        for line in xls_r.readlines():
            line = line.strip().split('\t')
            id = line[0]
            before = [int(x) for x in line[1:]]
            after = [str(x/(sum(before)+0.0001)) for x in before]
            per.write(id + '\t' + '\t'.join(after) + '\n')
    return

tree = lambda: defaultdict(tree)
first_atcg = tree()
loc_atcg = tree()
# first_atcg = dict(dict(dict()))
# loc_atcg = dict(dict(dict()))
orders=["A","G","C","T"]

with open(args.fasta,'r') as fasta:
    mirnas = fasta.read().lstrip('>').split('\n>')
    for mirna in mirnas:
        id = mirna.split('\n')[0]
        seq = ''.join(mirna.split('\n')[1:])
        name = id.split('_')[0].lower()
        time = id.split('_')[2].lstrip('x')
        if not first_atcg[name][len(seq)][seq[0]]:
            first_atcg[name][len(seq)][seq[0]] = 0
        first_atcg[name][len(seq)][seq[0]] += int(time)
        if len(seq) >= 18:
            for i in range(len(seq)):
                if not loc_atcg[name][i+1][seq[i]]:
                    loc_atcg[name][i + 1][seq[i]] = 0
                loc_atcg[name][i+1][seq[i]] += int(time)

for name in names:
    with open('%s_loc_bias.xls'%cfg.get('NAME', name), 'w') as loc_w, open('%s_first_bias.xls'%cfg.get('NAME', name), 'w') as length_w:
        loc_w.write('Location\tA\tG\tC\tU\n')
        length_w.write('Length(nt)\tA\tG\tC\tU\n')
        for i in range(int(args.min), int(args.max)+1):
            line = str(i) + '\t'
            for nt in orders:
                if first_atcg[name][i][nt]:
                    line = line + str(first_atcg[name][i][nt]) + '\t'
                else:
                    line = line + '0\t'
            line = line.strip('\t') + '\n'
            print(line)
            length_w.write(line)

        for i in range(1, int(args.max)+1):
            line = str(i) + '\t'
            for nt in orders:
                if loc_atcg[name][i][nt]:
                    line = line + str(loc_atcg[name][i][nt]) + '\t'
                else:
                    line = line + '0\t'
            line = line.strip('\t') + '\n'
            print(line)
            loc_w.write(line)
    change_to_per('%s_loc_bias.xls'%cfg.get('NAME', name))
    change_to_per('%s_first_bias.xls'%cfg.get('NAME', name))






