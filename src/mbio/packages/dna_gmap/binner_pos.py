# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.19

import re
import gzip
import argparse


def bin_var_pos(bin_info, pop_vcf, outfile):
    """
    bin marker对应的变异位点详情
    """
    bin_dict = {}
    with open(bin_info, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            item = line.strip().split("\t")
            chr = item[0].split("_")[0]
            if chr not in bin_dict.keys():
                bin_dict[chr] = {}
            print item[-1].split(";")
            bin_dict[chr][item[0]] = item[-1].split(";")
    chr_list = bin_dict.keys()
    if re.search(r".+.gz", pop_vcf):
        pop_vcf = gzip.GzipFile(pop_vcf, "r")
    else:
        pop_vcf = open(pop_vcf, "r")
    with open(outfile, "w") as w:
        w.write("Marker ID\tChr\tPos\tRef\tAlt\n")
        for line in pop_vcf:
            if line.startswith("#"):
                continue
            item = line.strip().split("\t")
            chr = item[0]
            pos = item[1]
            if chr in chr_list:
                for marker in bin_dict[chr].keys():
                    # print pos
                    # print bin_dict[chr][marker]
                    if pos in bin_dict[chr][marker]:
                        w.write(marker + "\t" + chr + "\t" + pos + "\t" + item[3] + "\t" + item[4] + "\n")
                        # bin_dict[chr][marker].pop(pos)
                        break


# bin_info = "/mnt/ilustre/users/sanger-dev/workspace/20180619/Single_bin_marker2/BinnerCalculate/bin_info_pos1.xls"
# pop_vcf = "/mnt/ilustre/users/sanger-dev/workspace/20180619/Single_bin_marker2/BinnerCalculate/pop.final.vcf"
# outfile = "/mnt/ilustre/users/sanger-dev/workspace/20180619/Single_bin_marker2/BinnerCalculate/bin_var.xls"
# bin_var_pos(bin_info, pop_vcf, outfile)

parser = argparse.ArgumentParser()
parser.add_argument("-bin", "--bin_marker", help="Total.bin.marker", required=True)
parser.add_argument("-vcf", "--pop_vcf", help="pop.final.vcf.gz", required=True)
parser.add_argument("-o", "--outfile", help="output file", required=True)
args = vars(parser.parse_args())


bin_var_pos(args["bin_marker"], args["pop_vcf"], args["outfile"])
