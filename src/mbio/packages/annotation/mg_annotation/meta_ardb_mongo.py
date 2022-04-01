# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
## fix by qingchen.zhang@20201009


import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config

class meta_ardb_anno(object):
    """
    宏基因ardb数据详细注释信息
    """
    def __init__(self):
        #self.client = MongoClient('mongodb://10.100.200.129:27017')
        #self.mongodb = self.client.sanger_biodb
        #self.mongodb = Config().mongo_client.sanger_biodb
        #self.ardb = self.mongodb.ardb
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.ardb = self.mongodb.ardb # 这里更新了resistance_group字段 20201009

    def ardb_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Query\tARG\tType\tIdentity(%)\tAlign_len\tAntibiotic_type\tRequirement\tClass\tClass_description\tAntibiotic_class\n')
            query_list = []
            done_list = []
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    #query = line[0]
                    query = ('_').join(line[5].split("_")[:-1])
                    #print line[0],"\n",query
                    ardb_ID = line[10]
                    #print ardb_ID
                    align_len = line[2]
                    act_iden = line[3]
                    detail = self.ardb.find_one({"arg": ardb_ID})
                    if detail:
                        iden = detail["idy"]
                        #print query,"\t",act_iden,"\t",iden
                        if float(act_iden) >= float(iden):
                            class_des = detail["class_des"]
                            resistance = detail["resistance"]
                            if isinstance(resistance,list) == True:
                                resistance = ';'.join(resistance)
                            type = detail["type"]
                            genebank = detail["arg"]
                            ardb_class = detail["class"]
                            resistance_group = detail["resistance_group"]
                            if isinstance(ardb_class,list) == True:
                                ardb_class = ';'.join(ardb_class)
                            require = detail["required_type"]
                            if  isinstance(require,list) ==True:
                                require = ';'.join(require)
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, genebank, type, act_iden, align_len, resistance, require, ardb_class, class_des, resistance_group))
                    else:
                        print "wrong ID:",ardb_ID


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o',metavar='[query_detail_table]',required=True,help='output file name')
    parser.add_argument('-iden',metavar='[identity]',help='input ardb identity cutoff')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_ardb_anno = meta_ardb_anno()
    meta_ardb_anno.ardb_by_mongo(align_table, anno_table)


