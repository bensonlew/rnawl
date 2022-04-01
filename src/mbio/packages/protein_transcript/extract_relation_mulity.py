#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/21 15:53
@file    : extract_relation_mulity.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""

import re
import copy
from multiprocessing import Pool
import sys
import time


def chunk_file(in_file, chunk_size=1000):
    line_list = []
    i = 0
    with open(in_file, 'rb') as fr:
        for num, seq in enumerate(fr):
            line_list.append(seq)
            i += 1
            if i == chunk_size:
                yield line_list
                i = 0
                line_list = []
        else:
            if line_list:
                yield line_list

def extract_relation(lines, protein_list):
    tmp_dict = dict()
    for line in lines:
        for pro in protein_list:
            if re.search('\t%s\t' % pro, line):
                line = line.strip().split('\t')
                gene = line[0]
                if gene not in tmp_dict:
                    tmp_dict[gene] = set()
                tmp_dict[gene].add(pro)
                break
    return tmp_dict

def multi_process(in_file, protein_list, p_num=10, chunk_size=100, ):
    g2p = dict()
    pool = Pool(p_num)
    results = list()
    for line_list in chunk_file(in_file, chunk_size):
        future = pool.apply_async(extract_relation, (line_list,protein_list))
        while len(results) > 1000:
            for fu in results:
                if fu.ready():
                    g2p.update(fu.get())
                    results.remove(fu)
        results.append(future)
    pool.close()
    pool.join()
    for fu in results:
        while not fu.ready():
            time.sleep(1)
        tmp = fu.get()
        for i in tmp:
            if i in g2p:
                g2p[i] = set(list(g2p[i]) + list(tmp[i]))
            else:
                g2p[i] = tmp[i]
        g2p.update(fu.get())
    return g2p


if __name__ == '__main__':
    pep_file = sys.argv[1]
    biomart_file = sys.argv[2]
    with open(pep_file, 'r') as pep_r:
        pep = pep_r.read().split('\n>')
        protein_list = [x.lstrip('>').split('\n')[0].strip().split(' ')[0].strip() for x in pep]
        for n in range(len(protein_list)):
            if u'.' in protein_list[n]:
                protein_list[n] = re.sub('\.\d+$', '', protein_list[n])
    tmp_list = copy.copy(protein_list)
    g2p = multi_process(biomart_file, tmp_list)
    with open('g2p.list', 'w') as list_w:
        genes = sorted(list(set(g2p.keys())))
        for gene in genes:
            str_i = gene + '\t'
            if gene in g2p:
                str_i += ';'.join(g2p[gene]) + '\t'
            else:
                str_i += '_' + '\t'
            list_w.write(str_i.strip('\t') + '\n')