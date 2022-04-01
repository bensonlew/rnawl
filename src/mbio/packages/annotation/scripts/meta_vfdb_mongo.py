# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'


import sys
# from pymongo import MongoClient
import argparse
from biocluster.config import Config


class meta_vfdb_anno(object):
    """
    meta宏基因vfdb数据详细注释信息
    "Gene\tType\tIdentity\tEvalue\tResistance\tRequirement\tClass\tClass_Description\n"
    """
    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.client = MongoClient('mongodb://192.168.10.187:27017')
        # self.client = MongoClient(Config().MONGO_BIO_URI)
        # self.mongodb = self.client.sanger_biodb
        # self.vfdb = self.mongodb.vfdb_new
        self.client = Config().get_mongo_client(mtype="metagenomic")
        self.db = client[Config().get_mongo_dbname("metagenomic")]
        self.vfdb = self.db.vfdb_new

    def vfdb_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Gene\tGi_number\tVfdbGene\tGene_description\tVFs\tVF_Function\tSpecies\tLevel1\tLevel2\n')
            query_list = []
            done_list = []
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    #query = line[0]
                    query = ('_').join(line[5].split("_")[:-1])
                    #print line[0],"\n",query
                    Subject_ID = line[10]
                    VFG_ID = Subject_ID.split("(")[0]
                    #print ardb_ID
                    evalue = line[1]
                    act_iden = line[3]
                    detail = self.vfdb.find_one({"VFG_ID":VFG_ID })
                    if detail:
                        Gi = detail["Gi_number"]
                        gene = detail["Gene"]
                        #if isinstance(resistance,list) == True:
                         #   resistance = ','.join(resistance)
                        gene_des = detail["Gene_des"]
                        species = detail["Origin"]
                        VF = detail["VFs"]
                        Function = detail["VF_Function"]
                        level1 = detail["Level1"]
                        level2 = detail["Level2"]
                        outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, Gi,gene, gene_des, VF,Function,species,level1,level2))
                    else:
                        print "wrong ID"


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o',metavar='[query_detail_table]',required=True,help='output file name')
    #parser.add_argument('-iden',metavar='[identity]',help='input ardb identity cutoff')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_vfdb_anno = meta_vfdb_anno()
    meta_vfdb_anno.vfdb_by_mongo(align_table, anno_table)

