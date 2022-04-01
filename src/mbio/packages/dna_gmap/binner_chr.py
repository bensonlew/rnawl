# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.15

import re
import gzip
import argparse


def bin_chr_stat(marker_file, pop_vcf, bin_marker, outfile):
    """
    得到marker文件里染色体或sca对应的pos
    """
    chr_dict, stat_dict, bin_dict = {}, {}, {}
    chr_list = []
    with open(bin_marker, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            item = line.strip().split("\t")
            marker_id = item[0]
            chr = marker_id.split("_")[0]
            if chr not in chr_list:
                bin_dict[chr] = 0
                chr_list.append(chr)
            bin_dict[chr] += 1
    with open(marker_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            marker_id = line.strip().split("\t")[0]
            marker_id1 = marker_id.split("_")
            chr = marker_id1[0]
            pos = marker_id1[1]
            if chr not in chr_dict.keys():
                chr_dict[chr] = []
                stat_dict[chr] = {"indel": 0, "snp": 0}
            chr_dict[chr].append(pos)
    if re.search(r".+.gz", pop_vcf):
        pop_vcf = gzip.GzipFile(pop_vcf, "r")
    else:
        pop_vcf = open(pop_vcf, "r")
    for line in pop_vcf:
        if line.startswith("#"):
            continue
        item = line.strip().split("\t")
        chr = item[0]
        pos = item[1]
        if chr in chr_list:
            if pos in chr_dict[chr]:
                ref = item[3]
                alt = item[4]
                if len(item[3]) == 1 and len(alt) == 1:
                    stat_dict[chr]["snp"] += 1
                else:
                    stat_dict[chr]["indel"] += 1
    pop_vcf.close()
    with open(outfile, "w") as w:
        w.write("#Chromosome ID\tBin Number\tSNP Number Per Bin\tIndel Number Per Bin\n")
        for chr in chr_list:
            bin_chr = bin_dict[chr]
            try:
                indel_chr = float(stat_dict[chr]["indel"]) / bin_chr
            except:
                indel_chr = 0
            try:
                snp_chr = float(stat_dict[chr]["snp"]) / bin_chr
            except:
                snp_chr = 0
            w.write(chr + "\t" + str(bin_chr) + "\t" + str(snp_chr) + "\t" + str(indel_chr) + "\n")


# marker_file = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/cp_binmaker/pop.filtered.marker"
# pop_vcf = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/cp_binmaker/pop.final.vcf.gz"
# bin_marker = "/mnt/ilustre/users/sanger-dev/workspace/20180612/Single_binner_calculate2/BinnerCalculate/Total.bin.marker"
# outfile = "/mnt/ilustre/users/sanger-dev/workspace/20180612/Single_binner_calculate2/BinnerCalculate/test.xls"
# bin_chr_stat(marker_file, pop_vcf, bin_marker, outfile)


parser = argparse.ArgumentParser()
parser.add_argument("-bin", "--bin_marker", help="Total.bin.marker", required=True)
parser.add_argument("-vcf", "--pop_vcf", help="pop.final.vcf.gz", required=True)
parser.add_argument("-marker", "--marker_file", help="pop.filtered.marker", required=True)
parser.add_argument("-o", "--outfile", help="output file", required=True)
args = vars(parser.parse_args())

bin_chr_stat(args["marker_file"], args["pop_vcf"], args["bin_marker"], args["outfile"])
