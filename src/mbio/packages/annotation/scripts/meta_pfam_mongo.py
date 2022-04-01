# -*- coding: utf-8 -*-
# __author__ = "qingchen.zhang" 2018/10/11

import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config

class meta_pfam_anno(object):
    """
    宏基因组pfam数据注释详细信息
    """
    def __init__(self):
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic",ref=True)]
        self.pfam = self.mongodb.pfam

    def pfam_by_mongo(self, align_table, anno_table):
        """
        从mongo参考库pfam注释
        :param align_table: hmmscan软件计算结果
        :param anno_table: pfam注释结果
        :return:
        """
        with open(align_table, "r") as infile, open(anno_table, "wb") as outfile:
            outfile.write("#Query\tPfam_Accession\tPfam_ID\tPfam_description\tType\tCLAN_id\tAlign_len\tIdentity(%)\n")
            for line in infile:
                line = line.strip().split("\t")
                if not "#" in line[0]:
                    #query = ('_').join(line[2].split("_")[:-1])
                    query = line[2]
                    pfam_acc = line[0]
                    act_iden = line[9]
                    act_iden = act_iden*100
                    act_len = int(line[8]) - int(line[7]) + 1
                    detail = self.pfam.find_one({"pfam_accession": pfam_acc})
                    if detail:
                        pfam_id = detail["pfam_id"]
                        type = detail["type"]
                        desc = detail["pfam_description"]
                        clan_id = detail.get("clan_accession")
                        if clan_id:
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, pfam_acc, pfam_id, desc, type, clan_id, act_len, act_iden))
                        else:
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, pfam_acc, pfam_id, desc, type, '-', act_len, act_iden))
                    else:
                        print "wrong ID"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='Input align table')
    parser.add_argument('-o', required=True, help= "output file name")
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_pfam_anno =meta_pfam_anno()
    meta_pfam_anno.pfam_by_mongo(align_table, anno_table)



































