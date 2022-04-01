# -*- coding: utf-8 -*-
import os,argparse
import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config
import pandas as pd


class mg_kegg_level(object):
    """
    meta宏基因kegg pathway分级丰度
    """
    def __init__(self):
        ko_path = Config().SOFTWARE_DIR + "/database/KEGG/metag_database/kegg_v94.2/kegg_ko_v94.2.xls"
        info_col = ["ko_name", "level1", "level2", "level3"]
        self.ko_info = pd.read_csv(ko_path, sep='\t', index_col=0)[info_col]

    def kegg_level(self, anno_table, outdir):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        levelout = outdir + "/kegg_level_anno.xls"
        anno_all = outdir + "/gene_kegg_anno_all.xls"
        with open(anno_table, 'r') as infile, open(levelout, 'wb') as outfile, open(anno_all, "wb") as outfile2:
            outfile.write('#Query\tlevel1\tlevel2\tlevel3\n')
            head = infile.next().strip()
            outfile2.write(head + "\t" +"KEGG_Name\tLevel1\tLevel2\tLevel3\n")
            for line in infile:
                line = line.strip()
                line1 = line.split("\t")
                if not "#Query" in line1[0]:
                    outfile2.write(line + "\t")
                    query = line1[0]
                    KO = line1[2]
                    detail = self.kegg.find_one({"ko_id":KO })
                    if KO in self.ko_info.index:
                        ko_name, le1, le2, le3 = list(self.ko_info.loc[KO, ])
                        outfile2.write("{}\t{}\t{}\t{}\n".format(ko_name, le1, le2, le3))
                        le1 = le1.replace('; ', ';').split(';')
                        le2 = le2.replace('; ', ';').split(';')
                        le3 = le3.replace('; ', ';').split(';')
                        for i in range(len(le1)):
                            outfile.write("{}\t{}\t{}\t{}\n".format(query, le1, le2, le3))
                    else:
                        outfile2.write("-\t-\t-\t-" + "\n")
                        print "wrong ID", KO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='KEGG level ')
    parser.add_argument('-i',metavar='[KEGG.gene.annotate]',required=True,help='KEGG.gene.annotate.xls')
    parser.add_argument('-o',metavar='[outputdir]',required=True,help='the output dir')
    args = parser.parse_args()
    anno_table = args.i
    outdir = args.o
    kegg_level = mg_kegg_level()
    kegg_level.kegg_level(anno_table, outdir)
