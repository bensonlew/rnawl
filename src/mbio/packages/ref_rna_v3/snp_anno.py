# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# modified by shicaiping at 20190322
from __future__ import division
import re


def snp_anno(variant_function, exonic_variant_function, combine_vcf, snp_stat):
    exonic_info = {}
    function_info = {}
    first_write_line = ['CHROM', 'START', 'END', 'REF', 'ALT', 'QUAL', 'ANNO', 'GENE(in or nearby)', 'Depth', 'MUT_type', 'MUT_info']
    samples = []
    sample_depths=[]
    sample_genotype=[]

    with open(exonic_variant_function, "r") as ef:
        for line in ef:
            if re.match(r"#", line):
                continue
            else:
                line = line.split("\t")
                key = line[3] + " " + line[4] + " " + line[6] + " " + line[7]
                if key not in exonic_info:
                    item1 = line[1]
                    item2 = line[2]
                    exonic_info[key] = [item1, item2]
                else:
                    continue

    with open(variant_function, "r") as vf:
        for n, line in enumerate(vf):
            if re.match(r"#", line):
                continue
            else:
                line = line.split()
                key = line[2] + " " + line[3] + " " + line[5] + " " + line[6]
                anno = line[0]
                gene = line[1]
                start = line[3]
                end = line[4]
                if key not in function_info:
                    function_info[key] = [anno, gene, start, end]
                else:
                    continue

    with open(combine_vcf, "r") as vcf, open(snp_stat, "w") as w:
        for n, line in enumerate(vcf):
            if re.match(r"##", line):
                continue
            elif re.match(r"#", line):
                title_line = line.strip().split()
                samples = [i.split(".")[0] for i in title_line[9:]]
                sample_depths=[i+"_sampledepth" for i in samples]
                sample_genotype=[i+"genotype" for i in samples]
                w.write("\t".join(first_write_line)+"\t")
                w.write("\t".join(samples)+"\t")
                w.write("\t".join(sample_genotype) +"\t"+"type"+"\n")
            else:
                line = line.split()
                total_dp = re.search("DP=[\d]*", line[7])
                total_dp = total_dp.group().split("=")[-1]
                # if line[7].startswith("INDEL"):
                if "INDEL" in line[7]:
                    type="indel"
                else:
                    type="snp"
                #mapping_quality = re.search("MQ=[\d]*", line[7])
                #mapping_quality = mapping_quality.group().split("=")[-1]
                qual = line[5]
                ## 增加SNP过滤条件，dp需大于等于样本个数*2，qual需大于等于20
                #if total_dp <= len(samples)*2 or mapping_quality <= 20:
                #if total_dp <= len(samples) * 2:
                #    continue
                chr = line[0]
                pos = line[1]
                ref = line[3]
                alt = line[4]
                all_alt = alt.split(",")
                i = 1
                mut_type=""
                mut_info=""
                start=""
                end=""
                gene=""
                anno=""
                gt=[]
                gt.append(ref)
                gt.extend(all_alt)
                sdepth = []
                sgenogype=[]
                for alt_allele in all_alt:
                    ln_mark = chr + " " + pos + " " + ref + " " + alt_allele
                    if ln_mark in exonic_info:
                        mut_type1 = exonic_info[ln_mark][0]
                        mut_info1 = exonic_info[ln_mark][1]
                    else:
                        mut_type1 = "."
                        mut_info1 = "."
                    mut_type=mut_type+mut_type1+"|||"
                    mut_info=mut_info+mut_info1+"|||"
                    if ln_mark in function_info:
                        anno1 = function_info[ln_mark][0]
                        gene1 = function_info[ln_mark][1]
                        #start = function_info[ln_mark][2]
                        #end = function_info[ln_mark][3]
                    else:
                        anno1 = "."
                        gene1 = "."
                    if i ==1:
                        start = line[1]
                        end = line[1]
                        gene=gene1
                        anno=anno1
                    i += 1
                mut_type=mut_type.strip("|||")
                mut_info=mut_info.strip("|||")
                for l in line[9:]:
                    infos = l.split(":")
                    gt_info=infos[0]
                    if gt_info =="./.":
                        gt_sample="./."
                    else:
                        gt1=int(gt_info.split("/")[0])
                        gt2 = int(gt_info.split("/")[1])
                        gt_sample=gt[gt1]+"/"+gt[gt2]
                    sgenogype.append(gt_sample)
                    if len(infos) == 5 or len(infos) == 3 or len(infos) == 7:
                        dp=infos[1]
                        sdepth.append(dp)
                    elif len(infos) == 6:
                        dp = infos[3]
                        sdepth.append(dp)
                    else:
                        sdepth.append("-")

                w.write("\t".join([chr, start, end, ref, alt, qual, anno, gene, total_dp, mut_type, mut_info]) + "\t")
                w.write("\t".join(sdepth)+ "\t"+"\t".join(sgenogype)+"\t"+type+ "\n")