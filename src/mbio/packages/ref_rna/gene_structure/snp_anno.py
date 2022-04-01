# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
import re


def snp_anno(variant_function, exonic_variant_function, snp_stat):
    anno_dict = {}
    reads_name_list = []
    first_write_line = ['CHROM', 'START', 'END', 'REF', 'ALT', 'READS_NUM', 'ANNO', 'GENE(in or nearby)', 'MUT_type', 'MUT_info']
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
                w.write("\t".join([line[2], line[3], line[4], line[5],line[6], line[-1], line[0], line[1], MUT_type, MUT_info]) + "\n")