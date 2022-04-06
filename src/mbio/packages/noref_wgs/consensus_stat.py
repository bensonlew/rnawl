## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "HONGDONG"
# last_modify:20190110

import re
import argparse
import os
import gzip
from collections import defaultdict


def consensus_stat(file, output, samplenum):
    cons_coverage = {}
    cons_samples = defaultdict(list)
    num = 0
    length = 0
    depth = 0
    coverage_10 = 0
    coverage_20 = 0
    coverage_50 = 0
    coverage_80 = 0
    with gzip.open(file, 'rb') as r:
        for line in r:
            if re.match('#.*', line):
                continue
            else:
                temp = line.strip().split("\t")
                num += 1
                length += len(temp[5])
                depth += len(temp[4].split(','))
                sample_num = []
                for sample in temp[4].split(','):
                    if sample.split('_')[0] not in sample_num:
                        sample_num.append(sample.split('_')[0])
                cons_coverage[temp[1]] = len(sample_num)
                cons_samples[temp[1]].append(len(sample_num) / float(samplenum))
    for key in cons_samples.keys():
        if cons_samples[key][0] >= 0.1:
            coverage_10 += 1
        if cons_samples[key][0] >= 0.2:
            coverage_20 += 1
        if cons_samples[key][0] >= 0.5:
            coverage_50 += 1
        if cons_samples[key][0] >= 0.8:
            coverage_80 += 1
    with open(os.path.join(output, 'consensus_stat.txt'), 'w') as w:
        w.write("Consensus Number\t{}\n".format(coverage_10))
        # w.write("Average Length\t{}\n".format(length / num))
        # w.write("Average Depth\t{}\n".format(depth / num))
        w.write("Average Coverage\t{}%\n".format(round(float(depth) / num / samplenum, 4) * 100))
        w.write("Consensus Coverage (20%)\t{}\n".format(coverage_20))
        w.write("Consensus Coverage (50%)\t{}\n".format(coverage_50))
        w.write("Consensus Coverage (80%)\t{}\n".format(coverage_80))
    print "consensus_stat统计完成--将计算consensus_coverage"
    with open(os.path.join(output, 'consensus_coverage.txt'), 'w') as w:
        w.write("#consensus_coverage\tconsensus_number\n")
        nums_ = defaultdict(list)
        final_coverage = defaultdict(list)
        for key in cons_coverage.keys():
            nums_[cons_coverage[key]].append(key)
        for key in nums_:
            final_coverage[key] = len(nums_[key])
        tuple = sorted(final_coverage.items())
        if int(samplenum) <= 10:
            for m in tuple:
                rate = float(m[0]) / float(samplenum) * 100
                w.write("{}%\t{}\n".format(str(round(rate, 2)), m[1]))
        else:
            cov_0_10, cov_10_20, cov_20_30, cov_30_40, cov_40_50, cov_50_60, cov_60_70, cov_70_80, cov_80_90, cov_90_100\
                = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            for m in tuple:
                rate = float(m[0]) / float(samplenum) * 100
                if 0 < rate < 10:
                    cov_0_10 += int(m[1])
                elif 10 <= rate < 20:
                    cov_10_20 += int(m[1])
                elif 20 <= rate < 30:
                    cov_20_30 += int(m[1])
                elif 30 <= rate < 40:
                    cov_30_40 += int(m[1])
                elif 40 <= rate < 50:
                    cov_40_50 += int(m[1])
                elif 50 <= rate < 60:
                    cov_50_60 += int(m[1])
                elif 60 <= rate < 70:
                    cov_60_70 += int(m[1])
                elif 70 <= rate < 80:
                    cov_70_80 += int(m[1])
                elif 80 <= rate < 90:
                    cov_80_90 += int(m[1])
                elif 90 <= rate <= 100:
                    cov_90_100 += int(m[1])
            w.write("0%-10%\t{}\n".format(cov_0_10))
            w.write("10%-20%\t{}\n".format(cov_10_20))
            w.write("20%-30%\t{}\n".format(cov_20_30))
            w.write("30%-40%\t{}\n".format(cov_30_40))
            w.write("40%-50%\t{}\n".format(cov_40_50))
            w.write("50%-60%\t{}\n".format(cov_50_60))
            w.write("60%-70%\t{}\n".format(cov_60_70))
            w.write("70%-80%\t{}\n".format(cov_70_80))
            w.write("80%-90%\t{}\n".format(cov_80_90))
            w.write("90%-100%\t{}\n".format(cov_90_100))

if __name__ == '__main__':
    """
    一个consensus的coverage_20指的是consensus中覆盖到的样本个数 除以样本总个数 大于20%的记录

    """
    parser = argparse.ArgumentParser(description="用于统计consensus结果(python consensus_stat.py "
                                                 "-i catalog.tags.tsv.gz -o ./ -n 2)")
    parser.add_argument("-i", "--infile", type=str, help="input catalog tags.tsv.gz")
    parser.add_argument("-n", "--sample_num", type=int, help="input sample's num")
    parser.add_argument("-o", "--outfile", type=str, help="outfile path")
    args = parser.parse_args()
    consensus_stat(args.infile, args.outfile, args.sample_num)
