# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180129 by shaohua.yuan


import os
import gc
import pandas as pd
import argparse
import numpy as np
import collections


class profile_select(object):
    def __init__(self):
        self.group = collections.OrderedDict()

    def run_select(self, profile, outfile, group_method, select_column="GeneID",
                   gene_list=None, group_file=None, sam=None, total="T", top=0,
                   merge="merge", trans=False, filter=None, scale=None, before_s=None):
        self.merge = merge
        self.filter = filter
        self.scale = scale
        self.before_s = before_s
        self.total = total
        self.group_method = int(group_method)
        self.top = top
        if trans:
            profiletable = pd.read_csv(profile, sep='\t', header=None,dtype=str)
            profiletable.set_index(0, inplace=True)
            profiletable = [profiletable.T, ]
        else:
            profiletable = pd.read_csv(profile, sep='\t', chunksize=1000000)
        if gene_list:
            gene_list = pd.read_csv(gene_list, sep='\t',
                                    usecols=[select_column])
            gene_list = set(gene_list[select_column])
        if group_file:
            self.group = self.get_group(group_file, merge=merge)
            self.sam_list = self.group.keys()
        else:
            self.sam_list = sam.split(',')
        if scale:
            ori_file = before_s
            scale_file = outfile
        else:
            ori_file = outfile
        out = []
        if scale:
            out_scale = []
        out_header = True
        write_mode = 'w'
        for select in profiletable:
            if set(select_column.split(',')) - set(select.columns):
                select = select[[select.columns[0],] + self.sam_list]
            else:
                select = select[select_column.split(',') + self.sam_list]
            try:
                if gene_list:
                    select = select[select[select_column].isin(gene_list)]
                    # gene_list -= set(select[select_column])
                tmp = self.by_one(select, ori_file, out_header, write_mode)
                if isinstance(tmp, pd.DataFrame) and not tmp.empty:
                    out.append(tmp)
                    select = tmp
                if scale:
                    scale_data = self.scale_data(select, scale)
                    scale_data = self.by_one(scale_data, scale_file, out_header, write_mode, True)
                    if isinstance(tmp, pd.DataFrame) and not scale_data.empty:
                        out_scale.append(scale_data)
                out_header = None
                write_mode = 'a+'
            except Exception as e:
                raise Exception(e)
            gc.collect()
        if top > 0:
            out_df = pd.concat(out)
            out_df = out_df.sort_values(['Total'], ascending=False)[0:top]
            if scale:
                scale_df = pd.concat(out_scale)
                scale_df = scale_df.sort_values(['Total'], ascending=False)[0:top]
            if total != 'T':
                out_df.drop('Total', 1, inplace=True)
                scale and scale_df.drop('Total', 1, inplace=True)
            scale and scale_df.to_csv(scale_file, sep='\t', index=True)
            out_df.to_csv(ori_file, sep='\t', index=False)

    def by_one(self, df, outfile, out_header, write_mode, scale=False):
        # df_total = df[self.sam_list].sum(axis=1)
        df1 = df[self.sam_list].apply(pd.to_numeric)
        df_total = df1.sum(axis=1)
        df['Total'] = df_total
        if self.group_method != 0:
            df = self.group_deal(df, self.group_method, self.merge)
        if self.filter == 'True':
            df = df[df_total > 0]
        if self.top > 0:
            df = df.sort_values('Total', ascending=False).head(self.top)
            return df
        else:
            if self.total != 'T':
                df.drop('Total', 1)
            if scale:
                df.to_csv(outfile, sep='\t', index=True,
                          header=out_header, mode=write_mode)
            else:
                df.to_csv(outfile, sep='\t', index=False,
                          header=out_header, mode=write_mode)

    def get_group(self, group_file, merge="merge"):
        with open(group_file, "rb") as f:
            f.next()
            for line in f:
                line = line.strip().split("\t")
                sample = line[0]
                if merge == "merge":
                    group_name = "Group_" + line[1]
                else:
                    group_name = line[1]
                self.group[sample] = group_name
        return self.group

    def scale(self, df, method):
        df = df.set_index(df.cloumns[0]).drop('Total', 1)
        df = df[self.sam_list]
        mean = df.mean(axis=1)
        if method == 'Ctr':
            scale_data = df.T - mean
        else:
            std = df.std(axis=1)
            if method == 'Par':
                scale_data = (df.T - mean) / np.sqrt(std)
            else:
                scale_data = (df.T - mean) / std
        return scale_data.T.fillna(0)

    def group_deal(self, df, method, merge):
        if int(method) == 1:
            groups = df[self.sam_list].groupby(self.group, axis=1).sum()
        elif int(method) == 2:
            groups = df[self.sam_list].groupby(self.group, axis=1).mean()
        elif int(method) == 3:
            groups = df[self.sam_list].groupby(self.group, axis=1).median()
        if merge == "merge":
            select = df.merge(groups, left_index=True, right_index=True,
                              how='outer')
        else:
            select = groups
        return select


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
    top = int(args.top)
    merge = args.merge_group
    before_s = os.path.join(args.odir or '.', "select_before_scale.xls")

    select.run_select(profile, outfile, group_method, select_column=args.sc,
                      gene_list=gene_list, group_file=group,
                      sam=sam, total=total, top=top, merge=merge,
                      trans=args.trans, filter=args.filter,
                      scale=args.scale, before_s=before_s)
