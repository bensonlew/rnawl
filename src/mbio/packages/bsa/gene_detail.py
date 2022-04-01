# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.27

import re
import argparse
import gzip


def find_gene_detail(index_file, chrom, start, end, outfile):
    """
    基因详情页查找变异位点信息，根据chrom, start, end对index_file进行筛选
    index_file: 0X-0X的index-calc.result.index
    """
    with open(index_file, "r") as f, open(outfile, "w") as w:
        for line in f:
            if line.startswith("#"):
                w.write(line)
            if line.startswith(chrom):
                item = line.strip().split("\t")
                if int(item[1]) >= start and int(item[1]) <= end:
                    w.write(line)

def find_gene_variant_detail(index_file, gene_id, outfile):
    """
    基因详情页查找变异位点信息，根据gene_id对index_file进行筛选
    index_file: 0X-0X的index-calc.result.variant
    """
    with open(index_file, "r") as f, open(outfile, "w") as w:
        for line in f:
            if line.startswith("#"):
                w.write(line)
            item = line.strip().split("\t")
            anno = item[-5].split(";")
            gene_ids = []
            for s in anno:
                gene_ids.append(s.split("|")[-1])
            if gene_id in gene_ids:
                w.write(line)

def find_gene_seq(ref_fa, chrom, start, end, outfile):
    """
    基因详情页查找基因序列信息，根据chrom, start, end对ref.fa进行筛选
    """
    seqs = ""
    le = 0
    if re.search(r".+.gz", ref_fa):
        infile = gzip.GzipFile(ref_fa, "r")
        lines = infile.readlines()
    else:
        lines = open(ref_fa).readlines()
    with open(outfile, "w") as w:
        for i in range(len(lines)):
            if re.search(">{}.*".format(chrom), lines[i]):
                for j in range(i+1, len(lines)):
                    seq = lines[j].split("\n")[0]
                    le += len(seq)
                    if le >= start and le < end:
                        if len(seq[start:]) == 0:
                            seqs += seq
                        else:
                            seqs += seq[start:]
                    if le >= start and le >= end:
                        mi = start - (le - len(seq))
                        ma = end - (le - len(seq))
                        seqs += seq[mi:ma+1]
                        if le == end:
                            seq = lines[j+1].split("\n")[0]
                            seqs += seq[0]
                        w.write(seqs)
                        break
                    if re.search(">.+", lines[j]):
                        break
                break


parser = argparse.ArgumentParser()
parser.add_argument("-t", "--type", help="method, detail/seq", required=True)
parser.add_argument("-i", "--index_file", help="index-calc.result.variant, 0X-0X", required=True)
parser.add_argument("-g", "--gene_id", help="gene_id", required=False)
parser.add_argument("-chr", "--chrom", help="chromosome", required=False)
parser.add_argument("-start", "--start", help="gene location of start", required=False)
parser.add_argument("-end", "--end", help="gene location of end", required=False)
parser.add_argument("-o", "--outfile", help="output file", required=True)
args = vars(parser.parse_args())

# if args["type"] == "detail":
#     find_gene_detail(args["index_file"], args["chrom"], int(args["start"]), int(args["end"]), args["outfile"])
if args["type"] == "detail":
    find_gene_variant_detail(args["index_file"], args["gene_id"], args["outfile"])
elif args["type"] == "seq":
    find_gene_seq(args["index_file"], args["chrom"], int(args["start"]), int(args["end"]), args["outfile"])
else:
    raise Exception("type类型错误")
