#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/28 16:13
@file    : extract_utr3_utr5.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import sys
import os
import textwrap
from collections import OrderedDict


def extract_utr3_utr5(all, part):
    utr3 = ''
    utr5 = ''
    run = True
    i=0
    l = len(part)
    while run:
        utr3 = all[0:i]
        all_ = all[i:]
        score = 0
        for n,nt in enumerate(part):
            if score < -1:
                i += 1
                break
            if score >= l - 2:
                run = False
                break
            if nt == all_[n]:
                score += 1
            else:
                score -= 1
    utr5 = all[len(utr3)+len(part):]
    return utr3,utr5


os.chdir('O:\\Users\\yitong.feng\\Desktop\lili\\Araport11_genes.201606.cdna&CDS.fasta')
all_file = 'Araport11_genes.201606.cdna.fasta'
part_file = 'Araport11_genes.201606.cds.fasta'

all_dict = OrderedDict()
with open(all_file) as ar:
    blocks = ar.read().split('\n')
    for block in blocks:
        block = block.strip().lstrip('>').split('\n')
        id = block[0].split(' ')[0]
        seq = ''.join(block[1:])
        all_dict[id] = seq

part_dict = OrderedDict()
with open(part_file) as pr:
    blocks = pr.read().split('\n')
    for block in blocks:
        block = block.strip().lstrip('>').split('\n')
        id = block[0].split(' ')[0]
        seq = ''.join(block[1:])
        part_dict[id] = seq

max_ = max([len(part_dict[p]) for p in part_dict])

with open('Araport11_genes_utr3.fasta', 'w') as u3, open('Araport11_genes_utr5.fasta', 'w') as u5:
    for p in part_dict:
        utr3, utr5 = extract_utr3_utr5(all_dict[p], part_dict[p], max_-10)
        u3.write('>' + p + '\n')
        u3.write('\n'.join(textwrap.wrap(utr3,width=60)) + '\n')
        u5.write('>' + p + '\n')
        u5.write('\n'.join(textwrap.wrap(utr5, width=60)) + '\n')