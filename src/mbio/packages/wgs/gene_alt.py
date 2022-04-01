# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.07

import re
import argparse
import gzip


def filter_pos(chrom, start, end, line):
    """
    根据start、end进行筛选
    """
    w_line = []
    jump = False
    item = line.strip().split("\t")
    if start <= int(item[1]) <= end:
        w_line.append(item[0])
        w_line.append(item[1])
        w_line.append(item[3])
        w_line.append(item[4])
        jump = True
    return jump, w_line


def get_gene_alt(gene_id, chrom, start, end, vcf_file, out_file):
    """
    通过基因位置，对vcf文件进行筛选，选出gene_id对应的alt信息
    """
    jump, jump_ = False, True
    if re.search(r".+.gz", vcf_file):
        infile = gzip.GzipFile(vcf_file, "r")
        lines = infile.readlines()
    else:
        lines = open(vcf_file).readlines()
    with open(out_file, "w") as w:
        w.write("#CHROM\tPOS\tREF\tALT\n")
        for line in lines:
            if line.startswith("#"):
                pass
            else:
                if line.startswith(chrom):
                    jump, w_line = filter_pos(chrom, start, end, line)
                    jump_ = jump
                    if jump:
                        w.write("\t".join(w_line) + "\n")
                if jump and not jump_:
                    break


parser = argparse.ArgumentParser()
parser.add_argument("-gene", "--gene_id", help="gene id", required=True)
parser.add_argument("-chr", "--chrom", help="chrom", required=True)
parser.add_argument("-start", "--start", help="the start of the gene_id", required=True)
parser.add_argument("-end", "--end", help="the end of the gene_id", required=True)
parser.add_argument("-vcf", "--vcf_file", help="the file of the pop final vcf", required=True)
parser.add_argument("-out", "--out_file", help="output file", required=True)
args = vars(parser.parse_args())


get_gene_alt(args["gene_id"], args["chrom"], int(args["start"]), int(args["end"]), args["vcf_file"], args["out_file"])
