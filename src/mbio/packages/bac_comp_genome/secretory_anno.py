# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
#20190917

import sys
import os,re
from pymongo import MongoClient
import argparse
from biocluster.config import Config


class meta_secretion_anno(object):
    """
    secretion数据详细注释信息
    """
    def secretion_table(self, ref, align_table, anno_table):
        anno = {}
        with open (ref, "r") as f:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                des = lin[1]+"\t"+lin[2]
                anno[lin[0]] =des
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Query\tGene\tType\tKO\tdes\n')
            for line in infile:
                line = line.strip().split("\t")
                if re.search("ko03070",line[4]):
                    gene_id = line[0]
                    ko =line[2]
                    des = line[3]
                    de = anno[line[2]]
                    outfile.write('{}\t{}\t{}\t{}\n'.format(gene_id, de, ko, des))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', metavar='[secretion_table]', required=True, help='Input datbase of secretion function')
    parser.add_argument('-i', metavar='[kegg_table]', required=True, help='Input kegg table')
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    meta_secretion_anno = meta_secretion_anno()
    meta_secretion_anno.secretion_table(args.s, args.i, args.o)
