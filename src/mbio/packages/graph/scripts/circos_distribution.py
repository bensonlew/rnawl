# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from __future__ import division
# import math
# from numpy import std,mean
import re
# import time
# import os
# import glob
import argparse


def get_length(length_file):
    chrom_length_dict = {}
    with open(length_file, "r") as f:
        for line in f:
            line = line.strip().split()
            # print line
            # chrom_length_dict[line[1]] = int(line[0])
            chrom_length_dict[line[0]] = {}
    print chrom_length_dict
    return chrom_length_dict


def count_reads_distribution(sam_file, chr_dict, range_num, outfile):
    count = 0
    new_chr_range_dict = chr_dict
    # print new_chr_range_dict
    with open(sam_file, "r") as f:
        # print f.readline()
        for line in f:
            if re.match(r"@", line):
                continue
            else:
                count += 1
                line = line.strip().split("\t")
                chr_name = line[2]
                # print chr_name
                # print line[3]
                pos = int(int(line[3])/range_num)
                ran = str(pos*range_num) + "-" + str((pos+1)*range_num)
                # print ran
                if chr_name in new_chr_range_dict:
                    if ran in new_chr_range_dict[chr_name]:
                        # print "cccccounting"
                        new_chr_range_dict[chr_name][ran] += 1
                    else:
                        # print "nnnnewwwww"
                        new_chr_range_dict[chr_name][ran] = 1
    for chr_stst in new_chr_range_dict:
        # print chr_stst
        with open(outfile + "/{}".format(chr_stst), "w") as w:
            # chrs = new_chr_range_dict.keys()
            w.write("range\t{}\n".format(chr_stst))
            for freq in new_chr_range_dict[chr_stst]:
                w.write("{}\t{}\n".format(freq, new_chr_range_dict[chr_stst][freq]))

    with open(outfile + "/all", "w") as w:
        w.write("{}".format(new_chr_range_dict))
    # print new_chr_range_dict

def main():
    _i = ''
    _l = ''
    _g = ''
    _o = ''
    parse = argparse.ArgumentParser()
    parse.add_argument('-i','--input_sam',help='input sam file')
    parse.add_argument('-l','--chr_lenght',help='chromsome length info')
    parse.add_argument('-f','--freq',help='base length')
    parse.add_argument('-o','--output_dir',help='output dir')
    args = parse.parse_args()
    _i = args.input_sam
    _l = args.chr_lenght
    _f = int(args.freq)
    _o = args.output_dir
    chr_length = get_length(_l)
    count_reads_distribution(_i, chr_length, _f, _o)
    # venn_table(_i,_g,_o)

if __name__ == '__main__':
    main()


# if __name__ == '__main__':
#     chr_length = get_length("/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/test_files/ref_rna/mapping/human_chr")
#     # print chr_length
#     # range_dict = range_chrom_length(chr_length, 100000)
#     # print range_dict
#     count_reads_distribution("/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/test_files/ref_rna/mapping/human.sam", chr_length, 100000, "/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/test_files/ref_rna/mapping/chr_stat")

