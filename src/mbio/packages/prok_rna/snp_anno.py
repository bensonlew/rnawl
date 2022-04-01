# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from __future__ import division
import re


def snp_anno(variant_function, exonic_variant_function, combine_vcf, snp_stat):
    # anno_dict = {}
    vcf_anno_dict = {}
    first_write_line = ['CHROM', 'START', 'END', 'REF', 'ALT', 'QUAL', 'ANNO', 'GENE(in or nearby)', 'Depth', 'MUT_type', 'MUT_info']
    samples = []
    with open(combine_vcf, "r") as vcf:
        for n, line in enumerate(vcf):
            if re.match(r"##", line):
                continue
            elif re.match(r"#", line):
                title_line = line.strip().split()
                samples = [i.split(".")[0] for i in title_line[9:]]
            else:
                line = line.split()
                ln_mark = line[0] + line[1]
                total_dp = re.search("DP=[\d]*", line[7])
                total_dp = total_dp.group()
                freqs = []
                qual = line[5]
                for l in line[9:]:
                    infos = l.split(":")
                    if len(infos) == 5:
                        freq = infos[1].split(",")[-1] + "/" + infos[2]
                        freqs.append(freq)
                    elif len(infos) == 6:
                        freq = infos[3].split(",")[-1] + "/" + infos[2]
                        freqs.append(freq)
                    else:
                        freqs.append("./.")
                # freqs包含很多样本的内容,这个嵌套列表里面有很多单元素的列表,这里下面先提前填充一个空元素,要不然下面的49,
                # 50行会把qual的面目全非
                vcf_anno_dict[ln_mark] = [total_dp, freqs, qual, []]

    with open(exonic_variant_function, "r") as ef:
        for line in ef:
            if re.match(r"#", line):
                continue
            else:
                line = line.split("\t")
                match_str = line[3] + line[4]
                # 通过染色体位置和postion合起来的字符串为键来判断文件,此处对外显子注释的那个文件加上MUT_type和ENST编号,
                # 下面往文件写入这些信息
                if match_str in vcf_anno_dict:
                    vcf_anno_dict[match_str][-1] = [line[1], line[2]]

    with open(variant_function, "r") as vf, open(snp_stat, "w") as w:
        w.write("\t".join(first_write_line)+"\t")
        w.write("\t".join(samples)+"\n")
        for n, line in enumerate(vf):
            if re.match(r"#", line):
                continue
            else:
                MUT_type = "."
                MUT_info = "."
                line = line.split()
                ln_mark = line[2] + line[3]
                # print ln_mark
                if ln_mark in vcf_anno_dict:
                    infos = vcf_anno_dict[ln_mark]
                    qual_value = vcf_anno_dict[ln_mark][-2]
                    w.write("\t".join([line[2], line[3], line[4], line[5], line[6], qual_value, line[0], line[1], infos[0].split("=")[-1]]) + "\t")
                    if len(vcf_anno_dict[ln_mark]) > 2:
                        if vcf_anno_dict[ln_mark][-1]:
                            MUT_type = vcf_anno_dict[ln_mark][-1][0]
                            MUT_info = vcf_anno_dict[ln_mark][-1][1]
                    del vcf_anno_dict[ln_mark]
                    w.write("\t".join([MUT_type, MUT_info]) + "\t")
                    w.write("\t".join(infos[1]) + "\n")
