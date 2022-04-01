# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'


import sys
from pymongo import MongoClient
import argparse
from biocluster.config import Config

class meta_card_anno(object):
    """
    meta宏基因card数据详细注释信息
    ardb.parse.anno.xls
    "Gene\tType\tIdentity\tEvalue\tResistance\tRequirement\tClass\tClass_Description\n"
    """
    def __init__(self):
        #self.client = MongoClient('mongodb://10.100.200.129:27017')
        #self.mongodb = self.client.sanger_biodb
        #self.mongodb = Config().biodb_mongo_client.sanger_biodb
        #更改mongodb连接方式
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]

        self.card = self.mongodb.card

    def card_by_mongo(self, align_table, anno_table):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Query\tARO_accession\tARO_name\tARO_description\tARO_category\tClass\n')
            #query_list = []
           #done_list = []
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    #query = line[0]
                    query = ('_').join(line[5].split("_")[:-1])
                    model_seq = line[10]
                    evalue = line[1]
                    act_iden = line[3]
                    detail = self.card.find_one({"model_seq": model_seq})
                    if detail:
                        ARO_accession = detail["ARO_accession"]
                        ARO_name = detail["ARO_name"]
                        ARO_category = detail["ARO_category"]
                        ARO_description = detail["ARO_description"]
                        Class = detail["Class"]
                        pass_value = detail["PASS_EVALUE"]
                        if float(evalue) <= 1e-30 :
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, ARO_accession, ARO_name, ARO_category,ARO_description,Class))
                    else:
                        print "wrong ID"


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o',metavar='[query_detail_table]',required=True,help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_card_anno = meta_card_anno()
    meta_card_anno.card_by_mongo(align_table, anno_table)
