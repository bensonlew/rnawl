# -*- coding: utf-8 -*-
from Bio import SeqIO
import pandas as pd

def circ_finish(genome_fa, genome_struction, output_path):
    genome = SeqIO.parse(genome_fa, 'fasta')
    genome_list = list()
    for i in genome:
        genome_list.append(i.name)
    try:
        with open(genome_struction, 'r') as g, open(output_path, 'w') as o:
            o.write('chr' + '\t' + 'start' + '\t' + 'end' + '\t' + 'direction' + '\t' + 'group' + '\t' + 'label' + '\n')
            for line in g.readlines()[1:]:
                chr, element_id, start, end, length = line.strip().split('\t')
                if int(start) <= int(end):
                    direction = 'p'
                    start_new = start
                    end_new = end
                else:
                    direction = 'n'
                    start_new = end
                    end_new = start
                o.write(chr + '\t' + str(start_new) + '\t' + str(end_new) + '\t' + direction + '\t' + "" + '\t' + element_id + '\n')
    except:
        g = pd.read_excel(genome_struction, header=0, sep='\t')
        with open(output_path, 'w') as o:
            o.write('chr' + '\t' + 'start' + '\t' + 'end' + '\t' + 'direction' + '\t' + 'group' + '\t' + 'label' + '\n')
            for i in g.index.tolist():
                chr = g.iloc[i]['Location']
                element_id = g.iloc[i]['Element ID']
                start = g.iloc[i]['Start']
                end = g.iloc[i]['End']
                length = g.iloc[i]['Length (bp)']
                if int(start) <= int(end):
                    direction = 'p'
                    start_new = start
                    end_new = end
                else:
                    direction = 'n'
                    start_new = end
                    end_new = start
                o.write(chr + '\t' + str(start_new) + '\t' + str(end_new) + '\t' + direction + '\t' + "" + '\t' + element_id + '\n')
def circ_scan(genome_fa, genome_struction, output_path):
    genome = SeqIO.parse(genome_fa, 'fasta')
    genome_dict = dict()
    for i in genome:
        genome_dict[i.name] = len(i.seq)
    a = sorted(genome_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    num = 0
    a_new = list()
    for j in a:
        a_new.append((j[0], num))
        num += j[1]
    a_new_dict = dict(a_new)
    try:
        with open(genome_struction, 'r') as g, open(output_path, 'w') as o:
            o.write('chr' + '\t' + 'start' + '\t' + 'end' + '\t' + 'direction' + '\t' + 'group' + '\t' + 'label' + '\n')
            for line in g.readlines()[1:]:
                chr, element_id, start, end, length = line.strip().split('\t')
                if int(start) <= int(end):
                    direction = 'p'
                    start_new = int(start) + int(a_new_dict[chr])
                    end_new = int(end) + int(a_new_dict[chr])
                else:
                    direction = 'n'
                    start_new = int(end) + int(a_new_dict[chr])
                    end_new = int(start) + int(a_new_dict[chr])
                o.write(chr + '\t' + str(start_new) + '\t' + str(end_new) + '\t' + direction + '\t' + "" + '\t' + element_id + '\n')
    except:
        g = pd.read_excel(genome_struction, header=0, sep='\t')
        with open(output_path, 'w') as o:
            o.write('chr' + '\t' + 'start' + '\t' + 'end' + '\t' + 'direction' + '\t' + 'group' + '\t' + 'label' + '\n')
            for i in g.index.tolist():
                chr = g.iloc[i]['Location']
                element_id = g.iloc[i]['Element ID']
                start = g.iloc[i]['Start']
                end = g.iloc[i]['End']
                length = g.iloc[i]['Length (bp)']
                if int(start) <= int(end):
                    direction = 'p'
                    start_new = int(start) + int(a_new_dict[chr])
                    end_new = int(end) + int(a_new_dict[chr])
                else:
                    direction = 'n'
                    start_new = int(end) + int(a_new_dict[chr])
                    end_new = int(start) + int(a_new_dict[chr])
                o.write(chr + '\t' + str(start_new) + '\t' + str(end_new) + '\t' + direction + '\t' + "" + '\t' + element_id + '\n')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="请输入基因组文件", type=str, required=True)
    parser.add_argument("-s", help="请输入基因组结构文件", type=str, required=True)
    parser.add_argument("-o", help="请输入输出文件路径", type=str, required=True)
    parser.add_argument("-t", help="请输入基因组类型", type=str, required=True)
    args = parser.parse_args()
    if args.t == 'finish':
        circ_finish(args.f, args.s, args.o)
    if args.t == 'scan':
        circ_scan(args.f, args.s, args.o)
