# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'


import sys
import os
from pymongo import MongoClient
import argparse
from biocluster.config import Config


class dna_card_anno(object):
    """
    dna card数据详细注释信息
    ardb.parse.anno.xls
    """

    def __init__(self):
        # self.client = MongoClient('mongodb://10.100.200.129:27017')
        # self.mongodb = self.client.sanger_biodb
        # self.mongodb = Config().biodb_mongo_client.sanger_biodb
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.card = self.mongodb.card

    def card_by_mongo(self, align_table, anno_table, category):
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write(
                'Gene ID\tARO_name\tARO_Accession\tARO_description\tARO_category\tEvalue\tIdentity(%)\tDrug_class\tResistance_mechanism\tScore\tCoverage(%)\n')
            # query_list = []
            # done_list = []
            self.aro = {}
            for line in infile:
                line = line.strip().split("\t")
                if not "Score" in line[0]:
                    # query = line[0]
                    query = line[5]
                    model_seq = line[10]
                    evalue = line[1]
                    align_len = line[2]
                    act_iden = line[3]
                    score = line[0]
                    coverge = round(abs(float(line[8])-float(line[7])) / float(line[6]), 3)
                    coverge = coverge * 100
                    if category== "aro_category":
                        detail = self.card.find_one({"model_seq": model_seq})
                    else:
                        detail = self.card.find_one({"aro_accession": model_seq})
                    if detail:
                        ARO_accession = detail["aro_accession"]
                        ARO_name = detail["aro_name"]
                        ARO_category = detail["aro_category"]
                        ARO_description = detail["aro_description"]
                        Class = detail["class"]
                        Class = Class.replace(",", ";")
                        class_des = detail["class_des"]
                        pass_value = detail["pass_evalue"]
                        drug_class = detail["drug_class"]
                        resistance_mechanism = detail["resistance_mechanism"]
                        if category== "aro_category":
                            category_name = ARO_category
                            head_name = "ARO Category"
                        else:
                            category_name = drug_class
                            head_name = "Drug class"
                        self.save_dict(self.aro, category_name)
                        outfile.write(
                            '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, ARO_name, ARO_accession,
                                                                                  ARO_description, ARO_category,
                                                                                  evalue, act_iden, drug_class,
                                                                                  resistance_mechanism, score, coverge))
                    else:
                        print
                        "wrong ID"
        level_file = os.path.dirname(anno_table) + '/card_category.xls'
        with open(level_file, 'wb') as outfile1:
            outfile1.write(head_name + '\tGene No.\n')
            for key in self.aro:
                if key != '-':
                    outfile1.write('{}\t{}\n'.format(key, self.aro[key]))

    def save_dict(self, mydict, category_name):
        for each in category_name.split(";"):
            if each in mydict:
                mydict[each] += 1
            else:
                mydict[each] = 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-c', metavar='[category]', help='aro_category or drug_class', default="aro_category")
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    category = args.c
    dna_card_anno = dna_card_anno()
    dna_card_anno.card_by_mongo(align_table, anno_table, category)