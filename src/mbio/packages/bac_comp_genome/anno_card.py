# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'


import sys
import os
import pandas as pd
import argparse
from biocluster.config import Config


class dna_card_anno(object):
    """
    dna card数据详细注释信息
    ardb.parse.anno.xls
    """
    def card_by_anno(self, ref ,align_table, anno_table):
        db = self.save_dict(ref)
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write(
                'Gene ID\tARO_name\tARO_Accession\tARO_description\tARO_category\tEvalue\tIdentity(%)\tDrug_class\tResistance_mechanism\tScore\tCoverage(%)\n')
            infile = infile.readlines()
            for line in infile[1:]:
                line = line.strip().split("\t")
                score = line[0]
                coverge = round(abs(float(line[8]) - float(line[7])) / float(line[6]), 3)
                coverge = coverge * 100
                query = line[5]
                evalue = line[1]
                act_iden = line[3]
                if line[10] in db.keys():
                    card_fun = db[line[10]].split("\t")
                    aro_name = card_fun[6]
                    aeo_accessipon = card_fun[0]
                    aro_des = card_fun[3]
                    aro_category = card_fun[5]
                    drug_class = card_fun[10]
                    resistance_mechanism =card_fun[11]
                    outfile.write(query + "\t" + aro_name + "\t" + aeo_accessipon+ "\t" + aro_des+ "\t" + aro_category+ "\t" + evalue+ "\t" + act_iden + "\t" + drug_class+ "\t" + resistance_mechanism+ "\t" + score+ "\t" + str(coverge) + "\n")
                else:
                    print 'line[10]'
                    print "wrong ID"

    def save_dict(self, file):
        dict = {}
        with open (file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                des = lin[0] + "\t" + "\t".join(lin[2:])
                dict[lin[1]] = des
        return dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]', required=True, help='Input xml table')
    parser.add_argument('-f', metavar='[function_database]', help='card function', default="function_database")
    parser.add_argument('-o', metavar='[query_detail_table]', required=True, help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    db = args.f
    dna_card_anno = dna_card_anno()
    dna_card_anno.card_by_anno(db, align_table, anno_table)