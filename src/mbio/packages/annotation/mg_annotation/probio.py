# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import sys
import argparse
import os
import pandas as pd

def probio_anno(annofile, ref_file, outfile):
    ref_dict = get_ref(ref_file)
    with open(annofile,"r") as f, open(outfile, "w") as outf:
        outf.write("#Query\tProbiotic_name\tBrand\tStrain\tGenus\tCommercial_Development_Stage\tUse_in\tProbiotic_Effect\tDisease_class\tICD_10_Disease_Code\tLineage\tProbio_ID\n")
        #outf.write("#GeneID\tProbio_ID\tName of probiotics\tStrain\n")
        for line in f:
            if not "#" in line:
                line = line.strip().split("\t")
                gene = line[0]
                tax_id = line[1]
                #genus = line[8]
                tax_list = [line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9]]
                lineage = ";".join(tax_list)
                if ref_dict.has_key(tax_id):
                    detail = ref_dict[tax_id]
                    detail = detail.split("\t")
                    probio_id = detail[0]
                    name = detail[1]
                    brand = detail[2]
                    strain = detail[3]
                    status =  detail[5]
                    use_in = detail[4]
                    effect = detail[6]
                    ICD = detail[7]
                    dis_class = detail[8]
                    genus= detail[-1]
                    # all = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(name,brand,strain,status,use_in,effect,ICD,lineage,probio_id)
                    all = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gene,name,brand,strain,genus,status,use_in,effect,dis_class,ICD,lineage,probio_id)
                    outf.write(all + "\n")

def get_ref(ref_file):
    ref_dict = {}
    with open(ref_file,"r") as f:
        for line in f:
            line = line.strip()
            line1 = line.split("\t")
            name = line1[-2]
            if name != "Not find":
                ref_dict[name] = line
    return ref_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[gene_nr_anno.xls]', required=True, help='Input gene nr annotation table')
    parser.add_argument('-ref', metavar='[probio reference file]', help='probio reference file')
    parser.add_argument('-o', metavar='[output file]', required=True, help='output file name')
    args = parser.parse_args()
    annofile = args.i
    ref_file = args.ref
    outfile = args.o
    probio_anno(annofile, ref_file,outfile)
