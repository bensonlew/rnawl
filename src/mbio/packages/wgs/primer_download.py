# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.06.07

import re
import argparse


def get_ssr_primer_download_file(infile, outfile):
    """
    得到ssr引物设计所需的下载文件，主要是将同一个SSR的引物放在一行，每对引物往后追加
    """
    with open(infile, "r") as f, open(outfile, "w") as w:
        head = f.readline().strip().split("\t")
        with open(infile, "r") as f1:
            first = f1.readline().strip().split("\t")
            for first_line in f1:
                first_line = first_line.strip().split("\t")
                try:
                    num = len(first_line[7].split(";"))
                except:
                    num = 0
                if num != 0:
                    break
            header = "#Chr\tSSR.nr\tSSR type\tSSR\tSize\tStart\tEnd"
            for i in range(num):
                header += "\tFORWARD PRIMER1 (5'-3')\tTm(℃)\tGC(%)\tLength(bp)\tREVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT size(bp)"
            w.write(header + "\n")
        for line in f:
            new_line = []
            item = line.strip().split("\t")
            w.write(item[0] + "\t" + item[1] + "\t" + item[2] + "\t" + item[3] + "\t" + item[4] + "\t" + item[5] + "\t" + item[6])
            try:
                for_primer = item[7].split(";")
                for_tm = item[8].split(";")
                for_gc = item[9].split(";")
                for_len = item[10].split(";")
                rev_primer = item[11].split(";")
                rev_tm = item[12].split(";")
                rev_gc = item[13].split(";")
                rev_len = item[14].split(";")
                product_size = item[15].split(";")
                for i in range(num):
                    w.write("\t" + for_primer[i] + "\t" + for_tm[i] + "\t" + for_gc[i] + "\t" + for_len[i] + "\t" + rev_primer[i]
                            + "\t" + rev_tm[i] + "\t" + rev_gc[i] + "\t" + rev_len[i] + "\t" + product_size[i])
            except:
                for i in range(num):
                    w.write("\t-\t-\t-\t-\t-\t-\t-\t-\t-")
            w.write("\n")


def get_primer_download_file(infile, outfile):
    """
    得到引物设计所需的下载文件，主要是将同一个SSR的引物放在一行，每对引物往后追加
    """
    with open(infile, "r") as f, open(outfile, "w") as w:
        head = f.readline().strip().split("\t")
        with open(infile, "r") as f1:
            first = f1.readline().strip().split("\t")
            for first_line in f1:
                first_line = first_line.strip().split("\t")
                try:
                    num = len(first_line[4].split(";"))
                except:
                    num = 0
                if num != 0:
                    break
            header = "#Chr\tPos\tRef\tAlt"
            for i in range(num):
                header += "\tFORWARD PRIMER1 (5'-3')\tTm(℃)\tGC(%)\tLength(bp)\tREVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT size(bp)\tVariation start(bp)\tVariation end(bp)"
            w.write(header + "\n")
        for line in f:
            new_line = []
            item = line.strip().split("\t")
            w.write(item[0] + "\t" + item[1] + "\t" + item[2] + "\t" + item[3])
            try:
                for_primer = item[4].split(";")
                for_tm = item[5].split(";")
                for_gc = item[6].split(";")
                for_len = item[7].split(";")
                rev_primer = item[8].split(";")
                rev_tm = item[9].split(";")
                rev_gc = item[10].split(";")
                rev_len = item[11].split(";")
                product_size = item[12].split(";")
                var_start = item[13].split(";")
                var_end = item[14].split(";")
                for i in range(num):
                    w.write("\t" + for_primer[i] + "\t" + for_tm[i] + "\t" + for_gc[i] + "\t" + for_len[i] + "\t" + rev_primer[i]
                            + "\t" + rev_tm[i] + "\t" + rev_gc[i] + "\t" + rev_len[i] + "\t" + product_size[i] + "\t" + var_start[i] + "\t" + var_end[i])
            except:
                for i in range(num):
                    w.write("\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-")
            w.write("\n")


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="infile", required=True)
parser.add_argument("-o", "--outfile", help="outfile", required=True)
parser.add_argument("-type", "--type", help="type", required=True)
args = vars(parser.parse_args())

if args["type"] not in ["ssr", "primer"]:
    raise Exception("type类型只能为ssr/primer")
if args["type"] == "ssr":
    get_ssr_primer_download_file(args["infile"], args["outfile"])
else:
    get_primer_download_file(args["infile"], args["outfile"])
