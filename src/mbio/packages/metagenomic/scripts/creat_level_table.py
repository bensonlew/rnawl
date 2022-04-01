# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import sys
import argparse
import os
import pandas as pd
from mbio.packages.statistical.large_data import table_parser


class CreateAbundTable(object):
    def __init__(self):
        self.anno_file = ""
        self.group = {}

    def select_table(self, anno_file, gene_list, profile, level, top, group_method, total, outfile, group_file=None,
                     database=None, sam=None, color_level=None,lowest_level=None, level_name=None):
        reader1 = pd.read_table(anno_file, sep='\t', header=0, iterator=True)
        anno_table = table_parser(reader1, chunk_size=100000, ignore_index=True).drop_duplicates('#Query')
        reader2 = pd.read_table(gene_list, sep='\t', header=0, iterator=True)
        gene_set = table_parser(reader2, chunk_size=100000, ignore_index=True)
        reader3 = pd.read_table(profile, sep='\t', header=0, iterator=True)
        profile_table = table_parser(reader3, chunk_size=100000, ignore_index=True).drop_duplicates('GeneID')
        # 根据gene_list筛选gene_anno 和profile
        profile_table = profile_table[profile_table["GeneID"].isin(list(gene_set["GeneID"]))]
        if sam != None:
            choose = "GeneID," + sam
            choose_list = choose.split(",")
            profile_table = profile_table.loc[:, choose_list]
            print profile_table.head()
        anno_table1 = anno_table[anno_table['#Query'].isin(list(gene_set["GeneID"]))]
        if len(profile_table) < 1:
            raise Exception('在所选profile集参数下数据为空，请重新设置该参数!')
        if len(anno_table) < 1:
            raise Exception('在所选gene_anno集参数下数据为空，请重新设置该参数!')
        #anno_table = anno_table1.ix[:, ["#Query", level]]
        #anno_table.columns = ["GeneID", level]
        if database == "nr":
            anno_table = anno_table1.ix[:, ["#Query"]]
            levels = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
            level_index = levels.index(level)
            anno_table[level] = anno_table1[levels[0]]
            if level_index > 0:
                for n in range(1, level_index + 1):
                    anno_table[level] += "," + anno_table1[levels[n]]
            anno_table.columns = ["GeneID", level]
        else:
            #anno_table = anno_table1.ix[:, ["#Query", level]]
            if level_name:
                anno_table = anno_table1.ix[:, ["#Query", lowest_level, level]]
                anno_table.columns = ["GeneID", lowest_level, level]
                color_level = None
            else:
                anno_table, color_level = self.get_color(database, anno_table1, level, color_level=color_level)
        print anno_table.head()
        abund = anno_table.merge(profile_table, on='GeneID', how='inner')
        abund[level]= abund[level].astype('str').replace("; ", ";")
        level_type_abund = abund.drop(level, axis=1).join(
            abund[level].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename(level))
        print level_type_abund.head()
        if color_level != None and database != "nr":
            level_type_abund[level] += "|" + level_type_abund[color_level]
            level_type_abund = level_type_abund.drop(color_level, axis=1)
        if level_name:
            level_type_abund_table = self.get_lowest(level_type_abund, lowest_level, level, level_name)
        else:
            level_type_abund = level_type_abund[(level_type_abund[level] != "-") & (level_type_abund[level] != "-|-")]
            level_type_abund_table = level_type_abund.groupby(level).sum()
        level_type_abund_table = level_type_abund_table.ix[list((level_type_abund_table > 0).any(axis=1))]
        print level_type_abund_table.head()
        if len(level_type_abund_table) < 1:
            raise Exception('在所选物种/功能的分类参数下数据为空，请重新设置该参数!')
        if group_file != None:
            level_type_abund_table = self.get_group_method(level_type_abund_table, group_method)
        if not "Total" in level_type_abund_table.columns:
            sample_list = level_type_abund_table.columns[1:len(level_type_abund_table.columns)]
            level_type_abund_table['Total'] = level_type_abund_table.loc[:, sample_list].apply(lambda x: x.sum(),
                                                                                               axis=1)
        level_type_abund_table = level_type_abund_table.sort_values(["Total"], ascending=False)
        #level_type_abund_table = level_type_abund_table.sort(["Total"], ascending=False)
        if top != "all":
            level_type_abund_table = self.get_top(gene_list, top, level_type_abund_table)
        if database == "nr":
            level_type_abund_table.index = level_type_abund_table.index.str.replace(",", ";")
        if int(total) == 0:
            level_type_abund_table = level_type_abund_table.drop(["Total"], axis=1)
        level_type_abund_table.to_csv(outfile, sep="\t")

    def get_lowest(self, abund, lowest_level, level_type, level_type_name):
        if lowest_level in ["Pathway"]:

            pathway = abund[lowest_level].str.replace("; ", ";").str.split(";", expand=True).stack().reset_index(level=1, drop=True).rename(lowest_level)
            level_type_abund = abund.drop(lowest_level, axis=1).join(pathway).reset_index(drop=True)
            level_type_abund = level_type_abund[level_type_abund[level_type] == level_type_name]
            level_type_abund = level_type_abund[level_type_abund[lowest_level] != "-"]
        else:
            level_type_abund = abund[abund[level_type] == level_type_name]
            level_type_abund = level_type_abund[level_type_abund[lowest_level] != "-"]
        level_type_abund_table = level_type_abund.groupby(lowest_level).sum()
        return level_type_abund_table

    def get_top(self, gene_list, top, level_type_abund_table):
        if len(gene_list) > int(top):
            level_type_abund_table = level_type_abund_table[0:int(top)]
        return level_type_abund_table

    def get_group_method(self,level_type_abund_table,group_method):
        self.group = self.get_group(group_file)
        print self.group
        level_type_abund_table = level_type_abund_table[self.group.keys()]
        if int(group_method) == 1:
            level_type_abund_table = level_type_abund_table.groupby(self.group, axis=1).sum()
        elif int(group_method) == 2:
            level_type_abund_table = level_type_abund_table.groupby(self.group, axis=1).mean()
        elif int(group_method) == 3:
            level_type_abund_table = level_type_abund_table.groupby(self.group, axis=1).median()
        else:
            level_type_abund_table = level_type_abund_table
        return level_type_abund_table

    def get_group(self, group_file):
        with open(group_file, "rb") as f:
            f.next()
            for line in f:
                line = line.strip().split("\t")
                sample = line[0]
                group_name = line[1]
                self.group[sample] = group_name
        return self.group

    def get_color(self, database, anno_ori, level, color_level=None):
        '''
        除NR以外的数据库颜色标注，同一个低水平对应两个高水平的话作为一个新颜色看待
        cog和kegg有一一对应关系，特殊处理，其他数据库直接带上选择的层级名称
        '''
        particular_list = ["card", "cazy", "ardb"]
        if color_level != None:
            if database == "cog" and color_level in ["Function","Category"]:
                anno_table = self.particular(anno_ori, level, color_level)
                color_level = None
            elif database == "kegg" and color_level in ["Level1","Level2"] :
                anno_table = self.particular(anno_ori, level, color_level)
                color_level = None
            else:
                anno_table = anno_ori.ix[:, ["#Query", level, color_level]]
                anno_table.columns = ["GeneID", level, color_level]
        else:
            anno_table = anno_ori.ix[:, ["#Query", level]]
            anno_table.columns = ["GeneID", level]
        return anno_table, color_level

    def particular(self, anno_ori, level, color_level):
        gene_t = anno_ori.ix[:, ["#Query"]]
        level_t = anno_ori[level].str.replace("; ", ";").str.split(';', expand=True).stack().reset_index(level=1, drop=True)
        color_level_t = anno_ori[color_level].str.replace("; ", ";").str.split(';', expand=True).stack().reset_index(level=1, drop=True)
        tmp_anno_table = pd.concat([level_t,color_level_t],axis=1)
        tmp_anno_table.columns = [level, color_level]
        tmp_anno_table2 = tmp_anno_table[level].str.cat(tmp_anno_table[color_level],sep="|")
        result = tmp_anno_table2.groupby(tmp_anno_table2.index).apply(';'.join)
        anno_table = pd.concat([gene_t, result],axis=1)
        anno_table.columns = ["GeneID", level]
        return anno_table



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[gene_*_anno.xls]', required=True, help='Input gene annotation table')
    parser.add_argument('-g', metavar='[geneset file]', help='Input gene set table')
    parser.add_argument('-p', metavar='[input profile]', help='Input profile table')
    parser.add_argument('-o', metavar='[output file]', required=True, help='output file name')
    parser.add_argument('-m', metavar='[map group]', help='group_file')
    parser.add_argument('-gm', metavar='[group sum method]', help='input group sum method', default=0)
    parser.add_argument('-l', metavar='level', help='input level')
    parser.add_argument('-database', metavar='database', help='input database')
    parser.add_argument('-t', metavar='top', help='input top', default="all")
    parser.add_argument('-cl', metavar='color_level', help='input color level')
    parser.add_argument('-ln', metavar='level name', help='input level name')
    parser.add_argument('-lowest', metavar='lowest level', help='input lowest level')
    #parser.add_argument('-database', metavar='[database]', required=True, help='input database name')
    parser.add_argument('-sam', metavar='[input selet columns]', help='input columns name,saun as samples')
    parser.add_argument('-total', metavar='[input is outfile need Total columns]', help='1 for True and 0 for False',
                        default=1)
    args = parser.parse_args()
    anno_file = args.i
    gene_list = args.g
    profile = args.p
    level = args.l
    outfile = args.o
    top = args.t
    group_method = args.gm
    total = args.total
    color_level = None
    database = None
    group_file = None
    sam = None
    level_name = None
    lowest_level = None
    if args.m:
        group_file = args.m
    if args.database:
        database = args.database
    if args.sam:
        sam = args.sam
    if args.cl:
        color_level = args.cl
    if args.ln:
        level_name = args.ln
        if not args.lowest:
            raise Exception('设置level_name时必须设置lowest level!')
        lowest_level = args.lowest
    run_select = CreateAbundTable()
    run_select.select_table(anno_file, gene_list, profile, level, top, group_method, total, outfile, group_file,
                            database, sam, color_level, lowest_level, level_name)
