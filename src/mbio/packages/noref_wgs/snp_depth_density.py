## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "HONGDONG"
# last_modify:20190110

import re
import argparse
import os
from collections import defaultdict


def tag_density(file, output):
    tag_num, num_tag, final_dict = defaultdict(list), defaultdict(list), defaultdict(list)
    sample_depth = defaultdict(list)
    sample_depth_sum = {}
    with open(file, 'rb') as r:
        for line in r:
            if re.match('##.*', line):
                continue
            elif re.match('#CHROM.*', line):
                item = line.strip().split("\t")
                sample_list = item[9:]
                for i in range(9, 9 + len(item[9:])):
                    sample_depth_sum[i] = {}
            else:
                item = line.strip().split("\t")
                if item[6] in ["PASS", "SNP", "INDEL", "FILTER"]:
                    tag_num[item[0]].append(1)  # tag \t snp num
                    for i in range(9, 9 + len(item[9:])):
                        try:
                            tem = int(item[i].split(":")[1])
                            if tem not in sample_depth[i]:
                                sample_depth_sum[i][tem] = 0
                                sample_depth[i].append(tem)
                            sample_depth_sum[i][tem] += 1
                        except:
                            pass
                else:
                    continue
    for key in tag_num.keys():
        num_tag[len(tag_num[key])].append(key)  # snp_num \t tag\
    for key in num_tag:
        final_dict[key] = len(num_tag[key])
    tuple = sorted(final_dict.items())
    with open(os.path.join(output, 'tag_snp.xls'), 'w') as w:
        w.write("#Tag_snp_num\ttag_number\n")
        for m in tuple:
            w.write("{}\t{}\n".format(m[0], m[1]))
    depth_list = []
    with open(os.path.join(output, 'depth_snp.xls'), 'w') as w:
        w.write("#Depth\t{}\n".format("\t".join(sample_list)))
        sum_ = []
        sample_max = {}
        for i in sample_depth.keys():
            depth_list = list(set(depth_list).union(set(sample_depth[i])))
            sum_.append(0)
            sample_max[i] = sorted(sample_depth[i])[-1]
        for dep in sorted(depth_list):
            for i in range(9, len(sample_list) + 9):
                if dep not in sample_depth_sum[i].keys():
                    if dep > sample_max[i]:
                        sum_[i-9] = None
                    else:
                        sum_[i-9] = sum_[i-9]
                else:
                    sum_[i-9] += sample_depth_sum[i][dep]
            w.write(str(dep) + "\t" + "\t".join([str(value) for value in sum_]) + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="用于统计tag中snp的个数")
    parser.add_argument("-i", "--infile", type=str, help="input populations.snps.vcf")
    parser.add_argument("-o", "--outfile", type=str, help="outfile path")
    args = parser.parse_args()
    tag_density(args.infile, args.outfile)
