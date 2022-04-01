# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190321

"""
合并ssr的结果，生成新的VCF
根据chr、pos相同做为同一个SSR进行合并
参考基因组：ref.fa.newmisa, chr1:1533	1	p1	(G)10	10	300	1543
样本SSR结果：ssr_result.list，每个样本对应的sample.ssr.result，#Chr    Pos     SSRbase Unit    Dep     Totaldep
vcf：#SSR ID  CHROM  Start  End  Repeat unit  Reference Repeat count  QUAL  INFO    FORMAT  sample1  sample2
     SSR ID为从1开始的整数编号;CHROM对应chr;Start对应pos;End若是在参考基因组上有这个SSR，则为start+size-1,没有为--;
     Repeat unit若是在参考基因组上有这个SSR，则为重复单元如AT，没有为--;Reference Repeat count若是在参考基因组上有这个SSR，则为size;
     QUAL为PASS，表示有样本和参考基因组都有这个SSR，为REF表示只有参考基因组有这个SSR，为ALT表示只有样本有这个SSR;
     INFO为参考基因组的SSR type;FORMAT为SB:ST:GT:AD:DP，对应样本里的SSRbase:SSRtype:Unit:Dep:Totaldep，需将样本里的冒号换成逗号；后面为每个样本的值
"""

import datetime
import argparse
from collections import defaultdict


def merge_ssr(ref_newmisa, ssr_list, vcf_path):
    chr_list = []
    ref_chr_pos, chr_pos = defaultdict(list), defaultdict(list)
    ref_dict = defaultdict(dict)
    with open(ref_newmisa, "rb") as f:  # chr1:1533	1	p1	(G)10	10	300	1543
        for line in f:
            if line.startswith("#") or line.startswith("ID"):
                continue
            item = line.strip().split("\t")
            chr = item[0].split(":")[0]
            pos = int(item[5])
            ref_chr_pos[chr].append(pos)
            ref_dict[chr][pos] = [item[2], item[3], item[4], item[5], item[6]]  # SSRtype SSRbase Unit Start End
    specimen_dict = {}
    specimen_list = []
    with open(ssr_list, "rb") as f:
        for line in f:
            item = line.strip().split("\t")
            sample = item[0]
            ssr_path = item[1]
            specimen_list.append(sample)
            sample_dict = defaultdict(dict)
            with open(ssr_path, "rb") as s:  # #Chr Pos SSRbase Unit Dep Totaldep
                for line in s:
                    if line.startswith("#") or line.startswith("ID"):
                        continue
                    item = line.strip().split("\t")
                    if int(item[-1]) < 10:
                        continue
                    unit_list = item[3].split(":")
                    dep_list = item[4].split(":")
                    unit_list_, dep_list_ = [], []
                    for i in range(len(unit_list)):
                        if int(unit_list[i]) * len(item[2]) < 8:
                            continue
                        if int(dep_list[i]) < 2:
                            continue
                        unit_list_.append(unit_list[i])
                        dep_list_.append(dep_list[i])
                    if not unit_list_:
                        continue
                    chr = item[0]
                    pos = int(item[1])
                    if chr not in chr_list:
                        chr_list.append(chr)
                    if pos not in chr_pos[chr]:
                        chr_pos[chr].append(pos)
                    sample_dict[chr][pos] = [item[2], ",".join(unit_list_), ",".join(dep_list_), item[5]]  # SSRbase Unit Dep Totaldep
            specimen_dict[sample] = sample_dict
    # print specimen_dict
    with open(vcf_path, "w") as w:
        w.write("#CHROM\tStart\tEnd\tSSR ID\tRepeat unit\tReference Repeat count\tQUAL\tINFO\tFORMAT\t{}\n".format("\t".join(specimen_list)))
        i = 0
        pass_num, ref_num, alt_num = 0, 0, 0
        for chr in chr_list:
            for pos in sorted(chr_pos[chr]):
                line = [chr, str(pos)]
                ref_ssr = False
                ref_info = ["--", str(i), "--", "--", "--"]
                if pos in ref_chr_pos:
                    ref_ssr = True
                    ref_info = [ref_dict[chr][pos][4], str(i), ref_dict[chr][pos][1], ref_dict[chr][pos][2], "ST=" + ref_dict[chr][pos][0]]
                line.extend(ref_info)
                line.append("SB:GT:AD:DP")
                sample_ssr = False
                for sample in specimen_list:
                    try:
                        sample_format = ":".join(specimen_dict[sample][chr][pos])
                        sample_ssr = True
                    except:
                        sample_format = "-:-:-:-"
                    line.append(sample_format)
                if ref_ssr and sample_ssr:
                    pass_num += 1
                    qual = "PASS"
                elif ref_ssr:
                    ref_num += 1
                    qual = "REF"
                elif sample_ssr:
                    alt_num += 1
                    qual = "ALT"
                else:
                    raise Exception(chr + str(pos) + "出错")
                line.insert(6, qual)
                i += 1
                if sample_ssr:
                    w.write("\t".join(line) + "\n")


parser = argparse.ArgumentParser(description="合并ref和各样本的SSR")
parser.add_argument("-r", "--ref_newmisa", help="ref.fa.newmisa", required=True)
parser.add_argument("-l", "--ssr_list", help="all sample ssr result list", required=True)
parser.add_argument("-o", "--vcf_path", help="output file", required=True)

args = vars(parser.parse_args())


merge_ssr(args["ref_newmisa"], args["ssr_list"], args["vcf_path"])
