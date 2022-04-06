## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "HONGDONG"
# last_modify:20190110

import re
import argparse
import os
import gzip
from collections import defaultdict


def tag_stat(file, output):
    tag_dict = defaultdict(list)
    tag_dep_dict = {}
    dep_tag = defaultdict(list)
    final_dict = {}
    with gzip.open(file, 'rb') as r:
        for line in r:
            if re.match('#.*', line):
                continue
            else:
                temp = line.strip().split("\t")
                if temp[2] not in ['consensus', 'model']:
                    tag_dict[temp[1]].append(1)
    for key in tag_dict.keys():
        tag_dep_dict[key] = len(tag_dict[key])
    for key in tag_dep_dict.keys():  # 每个深度中有多少个tag
        dep_tag[tag_dep_dict[key]].append(key)
    for key in dep_tag.keys():
        final_dict[key] = len(dep_tag[key])
    tuple = sorted(final_dict.items())
    with open(output + '_tag_dep.txt', 'w') as w:
        for m in tuple:
            w.write("{}\t{}\n".format(m[0], m[1]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="用于统计tag的深度分布")
    parser.add_argument("-i", "--infile", type=str, help="input sample's tags.tsv.gz")
    parser.add_argument("-o", "--outfile", type=str, help="outfile path")
    args = parser.parse_args()
    tag_stat(args.infile, args.outfile)
