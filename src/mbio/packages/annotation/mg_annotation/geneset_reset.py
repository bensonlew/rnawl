# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20171223


import sys
import argparse
import os
import pandas as pd
import shutil
import gc


class Geneset_reset(object):
    def __init__(self):
        self.sel_length_table = ""

    def run_select_profile(self, mytype, gene_length_table, profile, relative_profile, outDir, samples,
                           gene_file_list, select_gene_file= None, overview=None):
        self.remkdir(outDir)
        gene_dir = outDir + "/gene_profile"
        self.sel_length_table = gene_dir + "/gene.uniGeneset.fa.length.txt"
        uniGene_dir = outDir + "/uniGeneset"
        self.remkdir(gene_dir)
        self.remkdir(uniGene_dir)
        gene_list = []
        if select_gene_file:
            gene_table = pd.read_csv(select_gene_file, sep='\t', header=0, chunksize=5000000)
            for one in gene_table:
                gene_list.extend(one["GeneID"])
        if gene_file_list != "no_merge":
            files = gene_file_list.split(",")
            for eachfile in files:
                gene_table = pd.read_csv(eachfile, sep='\t', header=0, chunksize=5000000)
                for one in gene_table:
                    gene_list.extend(one["GeneID"])
        gene_list = set(gene_list)
        total_length, catalog_genes = self.table_select(gene_list, gene_length_table,
                                                        gene_dir + "/gene.uniGeneset.fa.length.txt",
                                                        stats=True)
        with open(uniGene_dir + "/geneCatalog_stat.xls", "w") as f:
            f.write("Catalog_genes\tCatalog_total_length(bp)\tCatalog_average_length(bp)\n")
            f.write("{}\t{}\t{}\n".format(catalog_genes, total_length,
                                          float(total_length) / catalog_genes))
        if samples != "all":
            choose_col = ["GeneID"] + samples.split(',')
        else:
            choose_col = None
        self.table_select(gene_list, profile, gene_dir + "/reads_number.xls",
                          choose_col=choose_col,
                          top_out=gene_dir + "/top100_reads_number.xls")
        self.table_select(gene_list, relative_profile, gene_dir + "/reads_number_relative.xls",
                          choose_col=choose_col,
                          top_out=gene_dir + "/top100_reads_number_relative.xls")
        if overview:
            if mytype == 2:
                outfile = outDir + "/anno_kegg.xls"
            else:
                outfile = outDir + "/anno_overview.xls"
            self.table_select(gene_list, overview, outfile, gene_col="#Query")

    def table_select(self, gene_list, file_name, outfile, gene_col="GeneID",
                     choose_col=None, top_out=None, stats=False):
        if os.path.exists(outfile):
            os.remove(outfile)
        df = pd.read_csv(file_name, sep='\t', header=0, chunksize=2000000, usecols=choose_col)
        tops = []
        total_length = 0
        catalog_genes = 0
        header = True
        for one in df:
            one = one[one[gene_col].isin(gene_list)]
            if choose_col:
                one["Total"] = one[choose_col[1:]].sum(axis=1)
            one.to_csv(outfile, sep='\t', index=False, mode='a+', header=header)
            header = False
            if stats:
                total_length += one["gene_length"].sum()
                catalog_genes += one["gene_length"].count()
            if top_out:
                one = one[one["Total"] != 0].sort_values("Total", ascending=False)[0:100]
                tops.append(one)
            del one
            gc.collect()
        if top_out:
            top_df = pd.concat(tops)[0:100]
            top_df.to_csv(top_out, sep='\t', index=False)
        if stats:
            return total_length, catalog_genes

    def run_length_distribute(self, ):
        len_dir = outDir + "/length_distribute"
        self.remkdir(len_dir)
        try:
            self.length_step(self.sel_length_table, 20, 200, len_dir + "/gene_step_200.final.txt")
            self.length_step(self.sel_length_table, 20, 300, len_dir + "/gene_step_300.final.txt")
            self.length_step(self.sel_length_table, 20, 400, len_dir + "/gene_step_400.final.txt")
            self.length_step(self.sel_length_table, 20, 500, len_dir + "/gene_step_500.final.txt")
            self.length_step(self.sel_length_table, 20, 600, len_dir + "/gene_step_600.final.txt")
            self.length_step(self.sel_length_table, 20, 800, len_dir + "/gene_step_800.final.txt")
        except Exception as e:
            raise Exception("创建基因长度分布表失败——{}".format(e))

    def remkdir(self, dir_name):
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)

    def length_step(self, gene_length_file, group_num, step, stat_out):
        with open(gene_length_file, "r") as r, open(stat_out, "a") as w:
            gene_length = []
            amount_group = []
            element_set = set("")
            for line in r:
                line = line.strip().split("\t")
                gene_len = line[1]
                if gene_len != "gene_length":
                    gene_length.append(gene_len)
            for f in gene_length:
                for i in range(group_num):
                    if (int(f) >= (i * step)) and (int(f) < ((i + 1) * step)):
                        amount_group.append(i)
                    else:
                        pass
                        # amount_group.append(group_num+1)
                    element_set.add(i)
            amount_group.sort()
            top_sum = 0
            for i in element_set:
                num_statistics = amount_group.count(i)
                if str(i) == '0':
                    area_line = str(i * step) + "~" + str((i + 1) * step) + "\t" + str(num_statistics) + "\n"
                    w.write(area_line)
                    top_sum += int(num_statistics)
                elif i < (group_num - 1):
                    area_line = str(i * step + 1) + "~" + str((i + 1) * step) + "\t" + str(num_statistics) + "\n"
                    w.write(area_line)
                    top_sum += int(num_statistics)
                else:
                    area_line = ">" + str(i * step) + "\t" + str(len(gene_length) - int(top_sum)) + "\n"
                    end_line = "total" + "\t" + str(len(gene_length)) + "\n"
                    w.write(area_line)
                    w.write(end_line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-gl', metavar='[gene_length_table]', required=True, help='Input gene_length_table')
    parser.add_argument('-p', metavar='[gene_profile]', required=True, help='Input gene_profile')
    parser.add_argument('-rp', metavar='[gene_relative_profile]', required=True, help='Input gene_relative_profile')
    parser.add_argument('-o', metavar='[output Directory]', required=True, help='output Directory name')
    parser.add_argument('-s', metavar='[select_genes_file]', help='Input select genes file')
    parser.add_argument('-sam', metavar='samples', help='input samples', default="all")
    parser.add_argument('-m', metavar='geneset_file_list', help='input geneset_file_list to merge with select genes',
                        default="no_merge")
    parser.add_argument('-overview', metavar='[overview file]', help='overview file')
    parser.add_argument('-ty', metavar='mytype', required=True, help='input type 2 for kegg')
    args = parser.parse_args()
    gene_length_table = args.gl
    samples = args.sam
    gene_file_list = args.m
    select_gene_file = None
    profile = args.p
    relative_profile = args.rp
    outDir = args.o
    mytype = args.ty
    overview = None
    if args.overview:
        overview = args.overview
    if args.s:
        select_gene_file = args.s
    if not args.s and not args.overview:
        raise Exception("不能同时没有select_gene_file和gene_list文件")
    run_reset = Geneset_reset()
    run_reset.run_select_profile(mytype, gene_length_table, profile, relative_profile, outDir, samples,
                                 gene_file_list, select_gene_file= select_gene_file, overview=overview)
    run_reset.run_length_distribute()
