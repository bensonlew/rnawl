#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/11 12:11
@file    : ucsc_refgene_extro_intro.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import os
import sys


def tr_perl(target_str, old, after):
    trans = ''
    for i in target_str:
        if i in old:
            i = after[old.index(i)]
        trans += i
    return trans

def extract_fa(seq, start, end, strand='+'):
    target = seq[start-1:end]
    if strand == '-':
        target = tr_perl(target, 'ATCGN', 'TAGCN')[::-1]
    return target

def extract_exon_intron_sites(ref_file, gene):
    starts = list()
    ends = list()
    chr_ = ''
    strand = '+'
    introns = list()
    with open(ref_file, 'r') as ref_r:
        for line in ref_r:
            if line.strip():
                line = line.strip().split('\t')
                gene_ = line[1]
                if gene_ == gene:
                    chr_ = line[2]
                    strand = line[3]
                    starts = line[9].strip(',').split(',')
                    ends = line[10].strip(',').split(',')
    if not starts:
        return
    starts = [int(x) for x in starts]
    ends = [int(x) for x in ends]
    exons = zip(starts,ends)
    for n, exon in enumerate(exons):
        try:
            s = exon[1] + 1
            e = exons[n+1][0] -1
            introns.append((s,e))
        except:
            pass
    return chr_, strand,exons, introns

if __name__ == '__main__':
    if len(sys.argv) != 5:
        exit('USAGE: %s ref_file fasta_dir gene out' % sys.argv[0])
    ref_file = sys.argv[1]
    fasta_dir = sys.argv[2]
    gene = sys.argv[3]
    out = sys.argv[4]
    chr_, strand, exons, introns = extract_exon_intron_sites(ref_file, gene)
    print(introns)
    seq = ''
    with open(os.path.join(fasta_dir, chr_+'.fa'), 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()
    with open('%s_exons.fasta'%out, 'w') as ew:
        for exon in exons:
            # print(exon[0], exon[1], strand)
            fa = extract_fa(seq, exon[0], exon[1], strand)
            ew.write('exon_%s_%s\n%s\n'%(str(exon[0]),str(exon[1]),fa))
    with open('%s_intron.fasta'%out, 'w') as iw:
        for intron in introns:
            fa = extract_fa(seq, intron[0], intron[1], strand)
            iw.write('intron_%s_%s\n%s\n'%(str(intron[0]),str(intron[1]),fa))




