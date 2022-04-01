# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last_modify:20190115
from sys import argv
import argparse
import os
import json
parse = argparse.ArgumentParser(description="snp_compare差异统计表计算")
parse.add_argument("-vcf", "--vcf_path", help="传入vcf文件", required=True)
parse.add_argument("-o", "--output", help="输出结果", required=True)
args = vars(parse.parse_args())
vcf_file = args["vcf_path"]
output_path = args["output"]

def stat_calc(vcf_file):
    with open(vcf_file) as f, open(output_path, "w") as w:
        lines=f.readlines()
        sum = 0  # 用于统计snp的个数
        num = 0  # 统计snp合法的个数
        sum_dep = 0  # 用于统计dep的求和
        miss = 0
        sample_num = 0
        for line in lines[1:]:
            if line.startswith("#"):
                pass
            else:
                sum = sum + 1
                row = line.strip().split("\t")
                sample_num =len(row[9:])
                for i in row[9:]:
                    try:
                        dep = i.strip().split(":")[1]
                        if dep != ".":
                            num = num + 1
                            sum_dep = sum_dep + int(dep)
                        else:
                            miss = miss + 1
                    except:
                        miss = miss + 1
        average_dep = float(sum_dep) / float(num)
        average_miss = float(miss) / float(int(sample_num) * int(sum))
        w.write("SNP Number" + "\t" + "Average Depth" + "\t" + "Miss Ratio" + "\n")
        w.write(str(sum) + "\t" + str(average_dep) + "\t" + str(average_miss)+ "\n")
stat_calc(vcf_file)
