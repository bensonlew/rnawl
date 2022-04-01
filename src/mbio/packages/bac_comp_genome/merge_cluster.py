# -*- coding: utf-8 -*-
# __author__ = 'shichang.xie'
## moodify qingchen.zhang

import re
import pandas as pd
import time, sys
import argparse


def id2cluster(cluster_file):
    id2c = {}
    with open(cluster_file, 'r') as r:
        for l in r:
            new_l = re.sub(r'\w+\||,|\s-\s', '\t', l)
            f = new_l.split()
            id2c.update({g: f[0] for g in f[3:]})
    return id2c


def merge(anno_file, id2c):
    anno = pd.read_csv(anno_file, sep='\t')
    name_list = ["Gene ID", "COG ID", "COG Description", "KO ID", "KO Description"]
    anno = anno[name_list]
    print 'get clsuter id column'
    anno['Cluster_ID'] = anno['Gene ID'].agg(lambda x: id2c[x] if x in id2c else '-')
    print 'merge'
    g_anno = anno.groupby('Cluster_ID').agg(lambda x: ';'.join(x)).reset_index()
    g_anno = g_anno[g_anno['Cluster_ID'] != '-']
    print 'done'
    g_anno.to_csv('t.xls', sep='\t', index=False)


def merge2(anno_file, id2c,out):
    anno = pd.read_csv(anno_file, sep='\t')
    name_list = ["Gene ID", "COG ID", "COG Description", "KO ID", "KO Description"]
    print 'get clsuter id column'
    anno['Cluster_ID'] = anno['Gene ID'].agg(lambda x: id2c[x] if x in id2c else '-')
    anno = anno[anno['Cluster_ID'] != '-']
    print 'merge'
    anno1 = anno[['Cluster_ID', "Gene ID"]].groupby('Cluster_ID').agg(lambda x: ';'.join(x)).reset_index()
    anno2 = anno[['Cluster_ID', "COG ID", "COG Description"]].drop_duplicates()
    anno2 = anno2[anno2['COG ID'] != '-']
    anno2 = anno2.groupby('Cluster_ID').agg(lambda x: ';'.join(x)).reset_index()
    new_anno = pd.merge(anno1, anno2, on='Cluster_ID', how='outer')

    anno3 = anno[['Cluster_ID', "KO ID", "KO Description"]].drop_duplicates()
    anno3 = anno3[anno3["KO ID"] != '-']
    anno3 = anno3.groupby('Cluster_ID').agg(lambda x: ';'.join(x)).reset_index()
    new_anno = pd.merge(new_anno, anno3, on='Cluster_ID', how='outer').fillna('-')
    new_anno.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', metavar='[annotation file]', required=True, help='anno file')
    parser.add_argument('-i', metavar='[homolog cluster]', required=True, help='homolog cluster file')
    parser.add_argument('-o', metavar='[out file]', required=True, help='result fasta')
    args = parser.parse_args()
    anno = args.a
    table = args.i
    out = args.o
    cluster_file, anno_file = table, anno
    print 'start'
    t1 = time.time()
    id2c = id2cluster(cluster_file)
    merge2(anno_file, id2c, out)
    t2 = time.time()
    print 'time: {}s'.format(t2 - t1)