# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180129 by shaohua.yuan


import os
import pandas as pd
import sys
import argparse
import numpy as np
import collections
from mbio.packages.statistical.large_data import table_parser


class profile_select(object):
    def __init__(self):
        self.group = collections.OrderedDict()

    def run_select(self, profile, outfile, group_method, select_column="GeneID", gene_list=None, group_file=None,
                   sam=None, total="T", top=0, merge="merge", trans=False, filter=None):
        if trans:
            reader = pd.read_table(profile, sep='\t',index_col=0, header=0, iterator=True)
            profiletable = table_parser(reader, chunk_size=100000, ignore_index=False)
            idx_name = profiletable.index.name
            profiletable = profiletable.T.reset_index().rename(columns={'index': idx_name})
        else:
            reader1 = pd.read_table(profile, sep='\t', header=0, iterator=True)
            profiletable = table_parser(reader1, chunk_size=100000, ignore_index=True)
        if gene_list:
            reader2 = pd.read_table(gene_list, sep='\t', header=0, iterator=True)
            gene_list = table_parser(reader2, chunk_size=100000, ignore_index=True)
            genes = gene_list[select_column]
            profiletable.index = profiletable[select_column]
            select = profiletable[profiletable[select_column].isin(genes)]
        else:
            select = profiletable
        select = select.drop_duplicates()
        if len(select) > 0:
            try:
                if group_file:
                    self.group = self.get_group(group_file, merge=merge)
                    if int(group_method) != 0:
                        if int(group_method) == 1:
                            groups = select[self.group.keys()].groupby(self.group, axis=1).sum()
                        elif int(group_method) == 2:
                            groups = select[self.group.keys()].groupby(self.group, axis=1).mean()
                        elif int(group_method) == 3:
                            groups = select[self.group.keys()].groupby(self.group, axis=1).median()
                        groups.index = select.index
                        ## add by shaohau.yuan 2018.3.29  group name is same with sample name
                        #groups_intersection = list(set(select.columns).intersection(set(groups.columns)))
                        #cols_to_use = list(set(groups.columns).difference(set(groups_intersection)))
                        #select = select.merge(groups[cols_to_use], left_index=True, right_index=True, how='outer')
                        groups.index = select.index
                        if merge == "merge":
                            select = select.merge(groups, left_index=True, right_index=True, how='outer')
                        else:
                            select = groups
                    else:
                        '''
                        if not trans:
                            tmp_index = select[select.columns[0]]
                        else:
                            tmp_index = select.index
                        '''
                        tmp_index = select[select.columns[0]]
                        select_group = []
                        [select_group.append(x) for x in self.group.keys() if x in select.columns]
                        select = select[select_group]
                        select.index = tmp_index
                    if int(top) > 0:
                        top = int(top)
                        sam1 = select.columns
                        select['Total'] = select.loc[:, sam1].apply(lambda x: x.sum(), axis=1)
                        select = select.sort_values(["Total"], ascending=False)[0:top]
                        if filter == 'True':
                            select = select[select["Total"] > 0]
                        select = select.drop('Total', 1)
                        #group_name = self.group.keys()
                        #select ['Total'] = select.loc[:, group_name].apply(lambda x: x.sum(), axis=1)
                if sam:
                    choose_names = select.columns[0] + "," + sam
                    choose_list = choose_names.split(",")
                    select = select.loc[:, choose_list]
                    if total == "T":
                        select['Total'] = select.loc[:, sam.split(",")].apply(lambda x: x.sum(), axis=1)
                        if filter == 'True':
                            select = select[select["Total"] > 0]
                        if int(top) > 0:
                            top = int(top)
                            select = select.sort_values(["Total"], ascending=False)[0:top]
                    else:
                        if int(top) > 0:
                            top = int(top)
                            select['Total'] = select.loc[:, sam.split(",")].apply(lambda x: x.sum(), axis=1)
                            select = select.sort_values(["Total"], ascending=False)[0:top]
                            if filter == 'True':
                                select = select[select["Total"] > 0]
                            select = select.drop('Total', 1)
                if merge != "merge":
                    select.to_csv(outfile, sep="\t", index=True)
                else:
                    select.to_csv(outfile, sep="\t", index=False)
            except Exception as e:
                raise Exception("选择{}失败！——{}".format(profile, e))
        else:
            select.to_csv(outfile, sep="\t", index=False)

    def get_group(self, group_file, merge="merge"):
        with open(group_file, "rb") as f:
            f.next()
            for line in f:
                line = line.strip().split("\t")
                sample = line[0]
                ## add by shaohau.yuan 2018.3.30  group name is same with sample name
                if merge == "merge":
                    group_name = "Group_" + line[1]
                else:
                    group_name = line[1]
                self.group[sample] = group_name
        return self.group

    def scale_data(self, ori_table, outDir, method='UV',sample=None):   #add method by zouguanqiing 20190605
        table = pd.read_table(ori_table, sep='\t', index_col=0)
        if sample:
            table = table[sample]   ##只用所选的样本做标准化 zouguanqing 2019
        if method=='Par':
            scaled_data = table.apply(lambda x: (x - np.mean(x)) / np.sqrt(np.std(x, ddof=1)), axis=1)
        elif method == 'Ctr':
            scaled_data = table.apply(lambda x: (x-np.mean(x)), axis=1)
        else:
            scaled_data = table.apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=1), axis=1)
        scaled_data =scaled_data.fillna(0)  # 防止处理后出现空值
        exp_profile = os.path.join(outDir, "scale_data.xls")
        scaled_data.to_csv(exp_profile, index=True, header=True, sep="\t")
        return exp_profile

    def orgin_abu_from_scale(self, origin_file, scale_top_file, outDir, group_method, group_file=None):
        table = pd.read_table(scale_top_file, sep='\t', index_col=0)
        select_names = table.index
        origin_table = pd.read_table(origin_file, sep='\t', index_col=0)
        selecl_origin = origin_table.loc[select_names,]
        select_abu = selecl_origin
        if not group_file:
            raise Exception("选择分组计算时没有分组文件")
        if not self.group:
            self.group = self.get_group(group_file, merge="nomerge")
        if int(group_method) == 0:
            select_group = []
            [select_group.append(x) for x in self.group.keys()]
            select_abu = select_abu[select_group]
        else:
            if int(group_method) == 1:
                groups = selecl_origin[self.group.keys()].groupby(self.group, axis=1).sum()
            elif int(group_method) == 2:
                groups = selecl_origin[self.group.keys()].groupby(self.group, axis=1).mean()
            elif int(group_method) == 3:
                groups = selecl_origin[self.group.keys()].groupby(self.group, axis=1).median()
            groups.index = selecl_origin.index
            select_abu = groups
        select_outfile = os.path.join(outDir, "select_before_scale.xls")
        select_abu.to_csv(select_outfile, sep="\t", index=True, header=True)


if __name__ == '__main__':
    # sam 方式挑选有Total，group_file方式挑选无Total
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[gene_profile]', required=True, help='Input gene profile')
    parser.add_argument('-s', metavar='[gene_selcet_list]', help='Input gene select_list')
    parser.add_argument('-o', metavar='[output file]', required=True, help='output file name')
    parser.add_argument('-sam', metavar='[select samples]', help='Input select samples')
    parser.add_argument('-g', metavar='[group file]',
                        help='input group file with two cloumns like samplename and group name!')
    parser.add_argument('-gm', metavar='[group method]', help='0 for 无，1 for sum, 2 for average,3 for middle!',
                        default=0)
    parser.add_argument('-st', metavar='[calculate total]', help='calculate total, T or F!', default="T")
    parser.add_argument('-top', metavar='[get top abundance]', help='get top abundance!', default=0)
    parser.add_argument('-sc', metavar='[select_column]', help='get select_column name!', default="GeneID")
    parser.add_argument('-merge_group', metavar='[merge_group]', help='get group or merge_group_file!', default="merge")
    parser.add_argument('-trans', metavar='[transposition table]', help='transposition table, T or F!', default=False)
    parser.add_argument('--scale', action='store_true', help="scale table")
    parser.add_argument('-odir', metavar='[outputDir]', help='outputDir for scale method')
    parser.add_argument('-filter', metavar='[filter]', help="whether filter 0")
    parser.add_argument('-scale_method', metavar='[scale_method]', help='scale method: UV, Ctr, Par', default='UV')
    args = parser.parse_args()
    profile = args.i
    outfile = args.o
    group_method = args.gm
    group = None
    sam = None
    gene_list = None
    if args.g:
        group = args.g
    if args.sam:
        sam = args.sam
    if args.s:
        gene_list = args.s
    select = profile_select()
    total = args.st
    top = args.top
    merge = args.merge_group
    if not args.scale:
        select.run_select(profile, outfile, group_method, select_column=args.sc, gene_list=gene_list, group_file=group,
                          sam=sam, total=total, top=top, merge=merge, trans=args.trans, filter=args.filter)
    else:
        '''
        先scale，后分组求和等，取top，生成scale和原丰度文件
        '''
        if group:
            group_data = pd.read_table(group,sep='\t',header=-1)
            sample_list = group_data[0].tolist()[1:]  #去掉表头
        elif sam:
            sample_list = sam.split(',')
        else:
            sample_list = None
        scale_file = select.scale_data(profile, args.odir,method=args.scale_method,sample=sample_list)  # 只用分析的样本做标准化
        select.run_select(scale_file, outfile, group_method, select_column=args.sc, gene_list=gene_list,
                          group_file=group, sam=sam, total=total, top=top, merge=merge, trans=args.trans)
        if group:
            select.orgin_abu_from_scale(profile, outfile, args.odir, group_method, group_file=group)
