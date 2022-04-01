# !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "yitong.feng"
# 20180721


import os
import sys
# from mako.template import Template
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

fna = sys.argv[1]
script_path = sys.argv[2]
bedtool_path = script_path + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'

operons = list()
transcripts = list()
rock_results = os.listdir('Rockhopper_Results')
seqs_op = set()
seqs_tr = set()

for file in rock_results:
    if u'transcripts' in file:
        transcripts.append(file)
    if u'operons' in file:
        seqs_op.add(file.split('_operons.txt')[0])
for file in transcripts:
    seqs_tr.add(file.split('_transcripts.txt')[0])
seqs_tmp = seqs_op & seqs_tr
seqs = set()

for i in os.listdir('rock_index'):
    if '.' in i:
        for j in seqs_tmp:
            if j in i:
                seqs.add(i)
    else:
        for j in seqs_tmp:
            if j == i:
                seqs.add(i)

seq_tr = set()


for file in transcripts:
    seqs_tr.add(file.split('_transcripts.txt')[0])


new_num = 0


with open('genome.predicted_RNA.bed', 'w') as genome_prebed:
    print("seqs is {}".format(seqs)) 
    for seq in seqs:
        with open('Rockhopper_Results/' + seq.split('.')[0] + '_transcripts.txt', 'r') as trans_r:
            header_trans = trans_r.readline()
            for line in trans_r.readlines():
                line = line.strip('\n').split('\t')
                if not len(line) < 6:
                    tran_start = line[0]
                    code_start = line[1]
                    code_end = line[2]
                    tran_end = line[3]
                    strand = line[4]
                    name = line[5]
                    type_ = line[6]
                    des = line[7]
                    if u'predicted' in type_:
                        if strand == '-':
                            tran_end, tran_start = tran_start, tran_end
                        length = int(tran_end) - int(tran_start) + 1
                        # id需要固定顺序以确保和blast结果一致
                        type__ = 'sRNA' + str(new_num).zfill(4)
                        new_num += 1
                        genome_prebed.write(
                            seq + '\t' + str(int(tran_start) - 1) + '\t' + tran_end + '\t' + type__ + '\t0' + '\t' + strand + '\t' + str(int(tran_start) - 1) + '\t' + tran_end + '\t0' + '\t1' + '\t' + str(length) + '\t0' + '\n')
#os.system('cat Rockhopper_Results/genome.gene.bed Rockhopper_Results/genome.knownnc.bed Rockhopper_Results/genome.predicted_RNA.bed >Rockhopper_Results/genome.feature.bed && cat Rockhopper_Results/fasta.gene.bed Rockhopper_Results/fasta.knownnc.bed Rockhopper_Results/fasta.predicted_RNA.bed >Rockhopper_Results/fasta.feature.bed')
# os.system('cat Rockhopper_Results/genome.gene.bed Rockhopper_Results/genome.predicted_RNA.bed >Rockhopper_Results/genome.feature.bed && cat Rockhopper_Results/fasta.gene.bed Rockhopper_Results/fasta.knownnc.bed Rockhopper_Results/fasta.predicted_RNA.bed >Rockhopper_Results/fasta.feature.bed')
bedtool_path = script_path + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'
cmd = '{} getfasta -fi {} -bed genome.predicted_RNA.bed -s -name -fo genome.predicted_RNA.fa'.format(
    bedtool_path, fna)
print cmd

os.system(cmd)

#  ------给genome.predicted_RNA.bed加入seqence那一列-----
