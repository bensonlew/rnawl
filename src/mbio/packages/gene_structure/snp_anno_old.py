# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from __future__ import division
import re


def snp_anno_test(variant_function, exonic_variant_function, combine_vcf, snp_stat):
    # anno_dict = {}
    vcf_anno_dict = {}
    first_write_line = ['CHROM', 'START', 'END', 'REF', 'ALT', 'ANNO', 'GENE(in or nearby)', 'Depth', 'MUT_type', 'MUT_info']
    samples = []
    with open(combine_vcf, "r") as vcf:
        for n, line in enumerate(vcf):
            if re.match(r"##", line):
                continue
            elif re.match(r"#", line):
                title_line = line.strip().split()
                samples = title_line[9:]
                print samples
            else:
                line = line.split()
                ln_mark =  line[0] + line[1]
                total_dp = re.search("DP=[\d]*", line[7])
                total_dp = total_dp.group()
                freqs = []
                for l in line[-4:]:
                    infos = l.split(":")
                    # print infos
                    if len(infos) == 5:
                        freq = infos[1].split(",")[-1] + "/" + infos[2]
                        # freq = infos[1].split(",")[-1] + "/" + infos[2]
                        freqs.append(freq)
                    else:
                        freqs.append("./.")
                vcf_anno_dict[ln_mark] = [total_dp, freqs]
    print len(vcf_anno_dict)

    with open(exonic_variant_function, "r") as ef:
        for line in ef:
            if re.match(r"#", line):
                continue
            else:
                line = line.split("\t")
                match_str = line[3] + line[4]
                if match_str in vcf_anno_dict:
                    vcf_anno_dict[match_str].append([line[1], line[2]])

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
                    # print infos[0]
                    w.write("\t".join([line[2], line[3], line[4], line[5], line[6], line[0], line[1], infos[0].split("=")[-1]]) + "\t")
                    # w.write("\t".join(infos[1]) + "\t")
                    # print "lllll"
                    # print vcf_anno_dict[ln_mark]
                    if len(vcf_anno_dict[ln_mark]) > 2:
                        MUT_type = vcf_anno_dict[ln_mark][-1][0]
                        # print MUT_type
                        MUT_info = vcf_anno_dict[ln_mark][-1][1]
                    del vcf_anno_dict[ln_mark]
                    w.write("\t".join([MUT_type, MUT_info]) + "\t")
                    w.write("\t".join(infos[1]) + "\n")


def snp_anno(variant_function, exonic_variant_function, snp_stat):
    anno_dict = {}
    reads_name_list = []
    first_write_line = ['CHROM', 'START', 'END', 'REF', 'ALT', 'ANNO', 'GENE(in or nearby)', 'MUT_type', 'MUT_info']
    with open(exonic_variant_function, "r") as ef:
        for line in ef:
            if re.match(r"#", line):
                continue
            else:
                line = line.split("\t")
                anno_dict[line[0]] = [line[1], line[2]]
    with open(variant_function, "r") as vf, open(snp_stat, "w") as w:
        w.write("\t".join(first_write_line)+"\n")
        for n, line in enumerate(vf):
            if re.match(r"#", line):
                continue
            else:
                line = line.split()
                ln_mark = "line" + str(n+1)
                if ln_mark in anno_dict:
                    MUT_type = anno_dict[ln_mark][0]
                    MUT_info = anno_dict[ln_mark][1]
                    del anno_dict[ln_mark]
                else:
                    MUT_type = "."
                    MUT_info = "."
                w.write("\t".join([line[2], line[3], line[4], line[5], line[6], line[0], line[1], MUT_type, MUT_info]) + "\n")


def snp_freq_stat(snp_vcf, annovar_stat, snp_out):
    freq_dict = {}
    first_line_list = ['CHROM', 'START', 'END', 'REF', 'ALT', 'READS_NUM', "MUT_rate", 'ANNO', 'GENE(in or nearby)', 'MUT_type', 'MUT_info']
    with open(snp_vcf, "r") as vcf:
        for line in vcf:
            if re.match(r"#", line):
                continue
            else:
                line = line.strip().split("\t")
                snp_key = ";".join([line[0], line[1], line[3], line[4]])
                if snp_key in freq_dict:
                    print "chongfuuuuuu"
                else:
                    values = line[-1].split(":")
                    if len(values) < 5:
                        continue
                    else:
                        freq_dict[snp_key] = int(values[1].split(",")[1])/int(values[2])
    # print len(freq_dict)
    with open(annovar_stat, "r") as f, open(snp_out, "w") as w:
        w.write("\t".join(first_line_list) + "\n")
        f.readline()
        for line in f:
            line = line.strip().split()
            snp_key = ";".join([line[0], line[1], line[3], line[4]])
            if snp_key in freq_dict:
                w.write("\t".join(line[:6]) + "\t" + str(freq_dict[snp_key]) + "\t" + "\t".join(line[6:10]) + "\n")
                del freq_dict[snp_key]


