# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'


import sys
# from pymongo import MongoClient
from biocluster.config import Config
import argparse

class meta_ardb_anno(object):
    """
    meta宏基因ardb数据详细注释信息
    ardb.parse.anno.xls
    "Gene\tType\tIdentity\tEvalue\tResistance\tRequirement\tClass\tClass_Description\n"
    """
    def __init__(self):
        #　self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.mongodb = self.client.sanger_biodb
        # self.ardb = self.mongodb.ardb
        self.client = Config().get_mongo_client(mtype="metagenomic")
        self.db = client[Config().get_mongo_dbname("metagenomic")]
        self.ardb = self.db.ardb


    def ardb_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            # anno.write('#Query\tCOG/NOG_Group\tFunction\tFunction_des\tFunction_type\tprotein_anno\tardb_des\n')
            outfile.write('#Gene\tGenBankID\tType\tIdentity\tEvalue\tResistance\tRequirement\tClass\tClassDescription\n')
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
                    evalue = line[1]
                    act_iden = line[3]
                    detail = self.ardb.find_one({"Subject_ID": ardb_ID})
                    if detail:
                        iden = detail["Idy"]
                        #print query,"\t",act_iden,"\t",iden
                        if float(act_iden) >= float(iden):
                            class_des = detail["Class_Info"]
                            resistance = detail["Resistance"]
                            if isinstance(resistance,list) == True:
                                resistance = ','.join(resistance)
                            type = detail["Type"]
                            genebank = detail["Subject_ID"]
                            ardb_class = detail["Class"]
                            if isinstance(ardb_class,list) == True:
                                ardb_class = ','.join(ardb_class)
                            require = detail["Required_type"]
                            if  isinstance(require,list) ==True:
                                require = ','.join(require)
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, genebank, type, act_iden,evalue,resistance,require,ardb_class,class_des))
                    else:
                        print "wrong ID"


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


