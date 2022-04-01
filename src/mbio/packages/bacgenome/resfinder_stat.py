# -*- coding: utf-8 -*-

import argparse
import os
from collections import defaultdict


def read_database(database):
    """
    将数据库的信息读进来
    :param database:
    phenotypes.txt
    :return:
    """
    gene_dict = defaultdict(list)
    with open(database, 'r') as f:
        f.readline()
        for line in f:
            spline = line.strip().split("\t")
            gene = spline[0].strip().split("_")[0]
            class_type = spline[1].strip()
            pheno = spline[2].strip()
            try:
                resistance = spline[4].strip()
            except:
                resistance = "-"
            all_anno = [gene, class_type, pheno, resistance]
            if gene not in gene_dict:
                gene_dict[gene] = [all_anno]
            else:
                gene_list = gene_dict[gene]
                if all_anno not in gene_list:
                    gene_list.append(all_anno)
                gene_dict[gene] = gene_list
    return gene_dict


def read_result(stat_file):
    """
    将resfinder的结果读进来并返回整理后的结果
    :param stat_file:
    :return:
    """
    gene_list = []
    anno_dict = defaultdict(list)
    with open(stat_file, 'r') as fr:
        fr.readline()
        for line in fr:
            spline = line.strip().split("\t")
            gene_iden = spline[1].strip()
            gene_cove = spline[3].strip()
            gene_name = spline[5].strip().split(" ")[0].strip()
            start = int(spline[5].strip().split(" ")[-2].strip())
            end = int(spline[5].strip().split(" ")[-1].strip())
            query_length = int(spline[2].strip().split("/")[0])
            hit_start = int(spline[4].strip().split("..")[0])
            hit_end = int(spline[4].strip().split("..")[1])
            query_start = int(spline[6].strip().split("..")[0])
            query_end = int(spline[6].strip().split("..")[1])
            if start > end:
                gene_start = end
                gene_end = start
            else:
                gene_start = start
                gene_end = end
            gene_start2 = gene_start + query_start
            gene_end2 = gene_start2 + query_length -1
            resis_gene = spline[0].strip()
            class_type = spline[7].strip()
            accession = spline[8].strip()

            all_anno = [gene_name, gene_start2, gene_end2, hit_start, hit_end, gene_iden, gene_cove, resis_gene, class_type, accession]
            if gene_name not in gene_list:
                gene_list.append(gene_name)
            if gene_name not in anno_dict:
                anno_dict[gene_name] = [all_anno]
            else:
                aa_list = anno_dict[gene_name]
                if all_anno not in aa_list:
                    aa_list.append(all_anno)
                anno_dict[gene_name] = aa_list
    return(gene_list, anno_dict)


def read_gff(gff):
    """
    将gff文件中的gene与location的对应关系确定
    :param gff:
    :return:
    """
    gff_dict = {}
    with open(gff, 'r') as fg:
        fg.readline()
        for line in fg:
            spline = line.strip().split("\t")
            gene = spline[0]
            location = spline[1].strip().split("_")[0]
            if gene not in gff_dict:
                gff_dict[gene] = location
    return gff_dict


def _main(input, output, database, sample, gff, type):
    """
    首先求出每个样本的最小序列长度
    然后根据最小序列长度过滤数据
    :param input: 输入文件
    :param output: 输出文件夹
    :param database : 需要挑选的database
    :param sample: 样品名称
    :param type: 数据库的类型，统计的class的类型方式不同
    :return:
    """
    if not os.path.exists(output):
        os.mkdir(output)
    stat_file = os.path.join(output, sample + ".stat.xls")
    class_file = os.path.join(output, sample+ "_" + type + ".class.xls")
    detail_file = os.path.join(output, sample + "_" + type + ".detail.xls")
    gene_dict = read_database(database)
    gene_list, anno_dict = read_result(input)
    gff_dict = read_gff(gff)
    with open(stat_file, 'w') as ws:
        ws.write("Sample Name\tGene No.\n")
        ws.write("{}\t{}\n".format(sample, str(len(gene_list))))
    with open(class_file, 'w') as wc, open(detail_file, 'w') as wd:
        wc.write("Sample Name\tClass\tGene No.\n")
        wd.write("Gene ID\tLocation\tSample Name\tResistance_gene\tClass\tPhenotype\tResistance Mechanism\tGene_Start\tGene_End\tHit_Start\tHit_End\tIdentity(%)\tCoverage(%)\tAccession no\n")
        class_dict = defaultdict(list)
        test_class_list = []
        for gene in gene_list:
            all_anno_list = anno_dict[gene]
            for all_anno in all_anno_list:
                gene_name = all_anno[0]
                gene_start = all_anno[1]
                gene_end = all_anno[2]
                hit_start = all_anno[3]
                hit_end = all_anno[4]
                gene_iden = all_anno[5]
                gene_cove = all_anno[6]
                resis_gene = all_anno[7]
                class_type2 = all_anno[8].split("resistance")[0].strip()
                accession = all_anno[9]
                if gene_name in gff_dict:
                    location = gff_dict[gene_name]
                else:
                    location = "-"
                if resis_gene in gene_dict:
                    all_class_list = gene_dict[resis_gene]
                    for class_list in all_class_list:##
                    # 这里判断逻辑是：同一个基因ID注释两个耐药基因，详情表就分两行展示，分类统计如果这两个耐药基因class一样就统计一个，如果class不一样，就分别统计，每个class一个
                    # 根据汪妍的需求修改
                        class_type = class_list[1]
                        print("aaaaa: {}".format(class_list))
                        if type in ['resfinder']:
                            class_type_list = class_list[1].split(",")
                            class_type_list = [x.strip() for x in class_type_list]
                            for c_type in class_type_list:
                                if c_type != "-" and c_type != "":
                                    if c_type not in class_dict:
                                        class_dict[c_type].append(gene_name)
                                    else:
                                        aa_list = class_dict[c_type]
                                        if gene_name not in aa_list:
                                            aa_list.append(gene_name)
                                        class_dict[c_type] = aa_list
                                if c_type not in test_class_list:
                                    test_class_list.append(c_type)
                        elif type in ['disinfinder']:
                            class_type_list = class_list[1].split(",")
                            class_type_list = [x.strip() for x in class_type_list]
                            for c_type in class_type_list:
                                if c_type != "-" and c_type != "":
                                    if c_type not in class_dict:
                                        class_dict[c_type].append(gene_name)
                                    else:
                                        aa_list = class_dict[c_type]
                                        if gene_name not in aa_list:
                                            aa_list.append(gene_name)
                                        class_dict[c_type] = aa_list
                                if c_type not in test_class_list:
                                    test_class_list.append(c_type)
                        pheno = class_list[2]
                        resistance = class_list[3]
                        if class_type != "-" and class_type != "":
                            aa = [gene_name, location, sample,resis_gene, class_type, pheno, resistance,gene_start, gene_end,hit_start, hit_end, gene_iden, gene_cove, accession]
                            bb = [str(x) for x in aa]
                            wd.write("{}\n".format("\t".join(bb)))
                        else:
                            aa = [gene_name, location, sample,resis_gene, class_type2, pheno, resistance,gene_start, gene_end,hit_start, hit_end, gene_iden, gene_cove, accession]
                            bb = [str(x) for x in aa]
                            wd.write("{}\n".format("\t".join(bb)))
        for cla in test_class_list:
            class_num = len(class_dict[cla])
            wc.write("{}\t{}\t{}\n".format(sample, cla, class_num))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[resfinder result]',required=True,help='resfinder result')
    parser.add_argument('-o', metavar='[output]', required=True, help='output_dir')
    parser.add_argument('-s', metavar='[database]', help='database')
    parser.add_argument('-n', metavar='[sample name]', help='sample name')
    parser.add_argument('-gff', metavar='[gff file]', help='Gff ile')
    parser.add_argument('-t', metavar='[type]', help='database type[resfinder, disinfinder]')
    args = parser.parse_args()
    input = args.i
    out_dir = args.o
    database = args.s
    sample = args.n
    gff = args.gff
    database_type = args.t
    _main(input, out_dir, database, sample, gff, database_type)

