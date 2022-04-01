# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from BCBio import GFF
import sys
import os
import re

gtf_file = sys.argv[1]
in_handle = open(gtf_file, 'r')
out_handle = open(gtf_file + ".name", 'w')

gene_id2name = dict()

'''
for rec in GFF.parse(in_handle):
    gene_id = rec.id
    gene_name = rec.name
    feature = rec.features
    types = [x.type for x in feature]
    if 'gene' in types:
        if gene_id2name.has_key(gene_id) and gene_id2name[gene_id] != gene_name:
            print "warning gene_id has diffrenent name {} {} {}".format(gene_id, gene_id2name[gene_id], gene_name)
        else:
            gene_id2name[gene_id] = gene_name

'''
for line in in_handle.readlines():
    cols = line.strip("\n").split("\t")
    # print cols[-1]
    tran_id = ""
    gene_name = ""
    m = re.match("transcript_id \"(.+)\";.*gene_name \"(.+)\";", cols[-1])
    if m:
        tran_id = m.group(1)
        gene_name = m.group(2)
    else:
        m = re.match("gene_name \"(.+)\";.*transcript_id \"(.+);\";", cols[-1])
        if m:
            tran_id = m.group(2)
            gene_name = m.group(1)

    if tran_id == "":
        continue

    if gene_id2name.has_key(tran_id) and gene_id2name[gene_name] != gene_name:
        print "warning gene_id has diffrenent name {} {} {}".format(tran_id, gene_id2name[tran_id], gene_name)
    else:
        gene_id2name[tran_id] = gene_name

for key, value in gene_id2name.items():
    out_handle.write("{}\t{}\n".format(key, value))
