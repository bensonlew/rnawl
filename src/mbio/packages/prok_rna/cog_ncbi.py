# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config
import time
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
from mbio.packages.rna.annot_config import AnnotConfig
import pandas as pd

class NcbiCogAnnot(object):
    """
    meta宏基因cog数据详细注释信息
    """
    def __init__(self, cog_version=2020):
        self.cog_version = cog_version
        self.import_db_file()

    def import_db_file(self):
        self.cog_files_dict = AnnotConfig().get_file_dict(db="cog", version=cog_version)
        seq2cog_table = self.cog_files_dict["cog-20.cog.csv"]
        self.seqid2cog = self.import_seqid2cog(seq2cog_table)

        cog_des_table = self.cog_files_dict["cog-20.def.tab"]
        self.cog2cat2, self.cog2des = self.import_cog_des(cog_des_table)

        cog_class_table= self.cog_files_dict["cog_description.xls"]
        self.cog_cat22des, self.cog_cat22cat1 = self.import_cog_class(cog_class_table)


    def import_seqid2cog(self, seq2cog_table):
        '''
        获取蛋白ID与cog对应关系
        '''
        seqid2cog = dict()
        with open(seq2cog_table, 'r') as f:
            for line in f:
                cols = line.strip().split(",")
                if cols[2] in seqid2cog:
                    seqid2cog[cols[2]].append([cols[6], cols[-1]])
                else:
                    seqid2cog[cols[2]] = [[cols[6], cols[-1]]]

        return seqid2cog

    def import_cog_class(self, cog_class_table):
        '''
        获取二级分类描述, 二级一级对饮关系
        '''
        a = pd.read_table(cog_class_table, sep="\t", header=0)
        cog_cat22des = dict(zip(a["Functional_categoryII_ID"], a['Functional_categoryII_description']))
        cog_cat22cat1 = dict(zip(a["Functional_categoryII_ID"], a['Functional_categoryI']))
        return cog_cat22des, cog_cat22cat1

    def import_cog_des(self, cog_des):
        '''
        获取cog 描述, cog与二级分类对应关系
        '''
        a = pd.read_table(cog_des, sep="\t", header=None)
        cog2cat2 = dict(zip(a[0], a[1]))
        cog2des = dict(zip(a[0], a[2]))
        return cog2cat2, cog2des


    def table2cog(self, align_table, anno_table):
        seq_cog_set = set()
        print("seqid2cog {}".format(self.seqid2cog.keys()[:5]))
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            infile.readline()
            outfile.write('#Query\tNOG\tNOG_description\tFunction\tFun_description\tCategory\tIdentity(%)\tAlign_len\n')
            for line in infile:
                cols = line.strip().split("\t")
                seq_id = cols[5]
                hit_id = cols[10]
                identity = cols[3]
                align_len = cols[2]
                
                if hit_id in self.seqid2cog:
                    # 仅判断是否比对上, 不考虑footprint
                    cogs = self.seqid2cog[hit_id]
                    # print cogs
                    for cog in cogs:
                        cog_id = cog[0]

                        if seq_id in seq_cog_set:
                            # 只保留一个cog
                            pass
                        else:
                            seq_cog_set.add(seq_id)
                            cog_des = self.cog2des[cog_id]
                            cog_cat2 = ";".join(self.cog2cat2[cog_id])
                            cog_cat2des = ";".join([self.cog_cat22des[cat] for cat in cog_cat2.split(";")])
                            cog_cat1 = ";".join([self.cog_cat22cat1[cat] for cat in cog_cat2.split(";")])
                            outfile.write("\t".join([
                                seq_id,
                                cog_id,
                                cog_des,
                                cog_cat2,
                                cog_cat2des,
                                cog_cat1,
                                identity,
                                align_len
                            ]) + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o',metavar='[query_detail_table]',required=True,help='output file name')
    parser.add_argument('-v',metavar='[cog version]',required=True,help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    cog_version = args.v
    cog_annot = NcbiCogAnnot(cog_version=cog_version)
    cog_annot.table2cog(align_table, anno_table)
