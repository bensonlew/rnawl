# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from optparse import OptionParser
import collections
import re
import os

parser = OptionParser(description='Get tabular file containing correspondence bewteen id, name and description')
parser.add_option('-i', '--input', dest='input', help='Input raw SNP annotation')
parser.add_option('-a', '--add', dest='add', help='The tabular file containing additional information')
parser.add_option('-o', '--output', dest='output', help='Output processed SNP annotation')
(opts, args) = parser.parse_args()

def main(snp_annotation_raw, id_name_des, snp_annotation_done):
    with open(snp_annotation_raw) as f1, open(id_name_des) as f2, open(snp_annotation_done, "w") as w:
        dict_id_symbol = collections.defaultdict(dict)
        dict_id_desc = collections.defaultdict(dict)
        _ = f2.readline()
        print 'INFO: start reading {}'.format(id_name_des)
        for line in f2:
            gene_id, symbol, des = line.strip().split("\t")
            dict_id_symbol[gene_id] = symbol
            dict_id_desc[gene_id] = des
        print 'INFO: start processing {}'.format(snp_annotation_raw)
        f1_col_names = f1.readline().strip().split("\t")
        w.write('{}\tGene name\tGene description\t{}\n'.format(
            '\t'.join(f1_col_names[0:8]), '\t'.join(f1_col_names[8:])
        ))
        for n, line in enumerate(f1):
            items = line.strip().split("\t")
            anno_list = items[6].split(";")
            geneby = items[7]
            gene_list = list()
            symbol_list = list()
            desc_list = list()
            if len(anno_list) >= 2:
                print 'INFO: find multiple ANNO of length {} at line {}'.format(len(anno_list), n + 1)
                if "splicing" in anno_list and "ncRNA_splicing" in anno_list and "intergenic" in anno_list and "ncRNA_intergenic" in anno_list:
                    print 'INFO: process line {} in case 01'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif "splicing" in anno_list and "ncRNA_splicing" in anno_list and "intergenic" in anno_list and ("ncRNA_intergenic" not in anno_list):
                    print 'INFO: process line {} in case 02'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                            else:
                                pass
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif "splicing" in anno_list and "ncRNA_splicing" in anno_list and ("intergenic" not in anno_list) and "ncRNA_intergenic" in anno_list:
                    print 'INFO: process line {} in case 03'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif "splicing" in anno_list and ("ncRNA_splicing" not in anno_list) and "intergenic" and "ncRNA_intergenic" in anno_list:
                    print 'INFO: process line {} in case 04'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                            else:
                                pass
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif ("splicing" not in anno_list) and "ncRNA_splicing" in anno_list and "intergenic" in anno_list and "ncRNA_intergenic" in anno_list:
                    print 'INFO: process line {} in case 05'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif "splicing" in anno_list and "ncRNA_splicing" in anno_list and ("intergenic" not in anno_list) and ("ncRNA_intergenic" not in anno_list):
                    print 'INFO: process line {} in case 06'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif "splicing" in anno_list and ("ncRNA_splicing" not in anno_list) and "intergenic" in anno_list and ("ncRNA_intergenic" not in anno_list):
                    print 'INFO: process line {} in case 07'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif "splicing" in anno_list and ("ncRNA_splicing" not in anno_list) and ("intergenic" not in anno_list) and "ncRNA_intergenic" in anno_list:
                    print 'INFO: process line {} in case 08'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif ("splicing" not in anno_list) and "ncRNA_splicing" in anno_list and "intergenic" in anno_list("ncRNA_intergenic" not in anno_list):
                    print 'INFO: process line {} in case 09'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif ("splicing" not in anno_list) and "ncRNA_splicing" in anno_list and ("intergenic" not in anno_list) and "ncRNA_intergenic" in anno_list:
                    print 'INFO: process line {} in case 10'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        elif anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif ("splicing" not in anno_list) and ("ncRNA_splicing" not in anno_list) and "intergenic" in anno_list and "ncRNA_intergenic" in anno_list:
                    print 'INFO: process line {} in case 11'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif ("splicing" not in anno_list) and ("ncRNA_splicing" not in anno_list) and ("intergenic" not in anno_list) and "ncRNA_intergenic" in anno_list:
                    print 'INFO: process line {} in case 12'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif ("splicing" not in anno_list) and ("ncRNA_splicing" not in anno_list) and "intergenic" in anno_list and ("ncRNA_intergenic" not in anno_list):
                    print 'INFO: process line {} in case 13'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "intergenic" or "ncRNA_intergenic":
                            tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(";")[i].split(",")]
                            tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                            gene_list.append(tmp_intergenic_list)
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif ("splicing" not in anno_list) and "ncRNA_splicing" in anno_list and ("intergenic" not in anno_list) and ("ncRNA_intergenic" not in anno_list):
                    print 'INFO: process line {} in case 14'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                elif "splicing" in anno_list and (not "ncRNA_splicing" in anno_list) and (not "intergenic" in anno_list) and (not "ncRNA_intergenic" in anno_list):
                    print 'INFO: process line {} in case 15'.format(n + 1)
                    for i in range(len(anno_list)):
                        if anno_list[i] == "splicing" or "ncRNA_splicing":
                            geneby_splicing = geneby.split(";")[i]
                            if ")," in geneby_splicing:
                                tmp_splicing_list_pre = [x.split("(") for x in geneby_splicing.split("),")]
                                tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                                gene_list.append(tmp_splicing_list)
                            elif ")," not in geneby_splicing and ")" not in geneby_splicing:
                                tmp_other_list = geneby_splicing.split(",")
                                gene_list.append(tmp_other_list)
                            elif ")," not in geneby_splicing and ")" in geneby_splicing:
                                for m in re.split(r"[(,]", geneby_splicing):
                                    if ":" not in m:
                                        gene_list.append([m])
                        else:
                            tmp_other_list = geneby.split(",")[i]
                            gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))

                else:
                    print 'INFO: process line {} in case 16'.format(n + 1)
                    tmp_other_list = re.split(r"[;,]", geneby)
                    gene_list.append(tmp_other_list)
                    for i in gene_list:
                        for j in i:
                            if j in dict_id_symbol.keys():
                                tmp_symbol = dict_id_symbol[j]
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = dict_id_desc[j]
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                            else:
                                tmp_symbol = "-"
                                symbol_list.append(tmp_symbol)
                                symbol_list_str = ";;".join(symbol_list)
                                tmp_desc = "-"
                                desc_list.append(tmp_desc)
                                desc_list_str = ";;".join(desc_list)
                    w.write('{}\t{}\t{}\t{}\n'.format(
                        '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                    ))
            else:
                if anno_list[0] == "splicing" or anno_list[0] == "ncRNA_splicing":
                    if ")," in geneby:
                        tmp_splicing_list_pre = [x.split("(") for x in geneby.split("),")]
                        tmp_splicing_list = [x[0] for x in tmp_splicing_list_pre]
                        gene_list.append(tmp_splicing_list)
                    elif ")," not in geneby and ")" not in geneby:
                        tmp_other_list = geneby.split(",")
                        gene_list.append(tmp_other_list)
                    elif ")," not in geneby and ")" in geneby:
                        for m in re.split(r"[(,]", geneby):
                            if ":" not in m:
                                gene_list.append([m])
                elif anno_list[0] == "intergenic" or anno_list[0] == "ncRNA_intergenic":
                    tmp_intergenic_list_pre = [x.split("(") for x in geneby.split(",")]
                    tmp_intergenic_list = [x[0] for x in tmp_intergenic_list_pre]
                    gene_list.append(tmp_intergenic_list)
                else:
                    tmp_other_list = geneby.split(",")
                    gene_list.append(tmp_other_list)
                symbol_list_str = "-"
                desc_list_str = "-"
                for i in gene_list:
                    for j in i:
                        if j in dict_id_symbol.keys():
                            tmp_symbol = dict_id_symbol[j]
                            symbol_list.append(tmp_symbol)
                            symbol_list_str = ";;".join(symbol_list)
                            tmp_desc = dict_id_desc[j]
                            desc_list.append(tmp_desc)
                            desc_list_str = ";;".join(desc_list)
                        else:
                            tmp_symbol = "-"
                            symbol_list.append(tmp_symbol)
                            symbol_list_str = ";;".join(symbol_list)
                            tmp_desc = "-"
                            desc_list.append(tmp_desc)
                            desc_list_str = ";;".join(desc_list)
                if len(gene_list) == 0:
                    print 'INFO: find empty gene list at line {}'.format(n + 1)
                w.write('{}\t{}\t{}\t{}\n'.format(
                    '\t'.join(items[0:8]), symbol_list_str, desc_list_str, '\t'.join(items[8:])
                ))
    if os.path.getsize(snp_annotation_done):
        print 'INFO: succeed in exporting {}'.format(snp_annotation_done)

if __name__ == '__main__':
    if opts.input and opts.add and opts.output:
        main(opts.input, opts.add, opts.output)
    else:
        parser.print_help()