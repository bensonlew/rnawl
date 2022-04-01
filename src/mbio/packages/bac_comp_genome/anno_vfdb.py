# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
#20190917

import sys
import os
from pymongo import MongoClient
import argparse
from biocluster.config import Config


class meta_vfdb_anno(object):
    """
    dna vfdb数据详细注释信息
    """
    def vfdb_table(self, ref, align_table, anno_table):
        anno = {}
        with open (ref, "r") as f:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split("\t")
                des = lin[1]+"\t"+lin[2]+"\t"+lin[3]+"\t"+lin[4]+"\t"+lin[5]+"\t"+lin[6]+"\t"+lin[7]+"\t"+lin[8]
                anno[lin[0]] =des
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Query\tVFDB ID\tGi_num\tVfdbGene\tGene_description\tSpecies\tVFs\tVFs_Description\tLevel2\tLevel1\tIdentity(%)\tEvalue\tScore\tCoverage(%)\n')
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    query = line[5]
                    Subject_ID = line[10]
                    align_len = line[2]
                    VFG_ID = Subject_ID.split("(")[0]
                    evalue = line[1]
                    act_iden = line[3]
                    coverge = round(abs(float(line[8])-float(line[7]))/float(line[6]), 3)
                    coverge = coverge*100
                    if VFG_ID in anno:
                        outfile.write(
                            '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, VFG_ID, anno[VFG_ID], act_iden, evalue, align_len, coverge))
                    else:
                        print
                        "wrong ID"

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', metavar='[anno_table]', required=True, help='Input datbase function')
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    meta_vfdb_anno = meta_vfdb_anno()
    meta_vfdb_anno.vfdb_table(args.s, args.i, args.o)
