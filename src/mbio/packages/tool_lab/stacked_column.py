# coding=utf-8
# __author__ = 'fwy'


from Bio import SeqIO
import fileinput
import re
import os
import subprocess
import urllib2
from collections import defaultdict
import regex
import pandas as pd
# from collections import Counter
# import matplotlib.pyplot as plt
import sys


def cat_samples_percent(exp,sep):
    sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[sep]
    table = pd.DataFrame(pd.read_table(exp, sep=sep, index_col=0))
    table.columns = [str(i) for i in table.columns]
    table_name = table.index.name
    table = table.loc[(table > 0).any(axis=1), :] #去除全为0的行
    abund = table
    abund['Col_sum'] = abund.apply(lambda x: x.sum(), axis=1)
    abund_table = abund.sort_values(by=['Col_sum'], ascending=0)
    del abund_table["Col_sum"]#这两部是根据总量进行排序
    abund_table_path = os.path.join("taxa.table.xls")
    abund_table.to_csv(abund_table_path, sep="\t", encoding="utf-8")
    abund_table.columns = [str(i) for i in abund_table.columns]
    abund_table.loc['Row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
    sample_empt = []
    b = abund_table.apply(lambda x: x.sum(), axis=0)
    for i in range(len(abund_table.columns)):
        if b[i] > 0:
            pass
        else:
            sample_name = abund_table.columns[i]
            sample_empt.append(sample_name)
    abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['Row_sum'], axis=1).drop('Row_sum')
    abund_table_percent.to_csv("taxa.percents.table.xls", sep="\t", encoding="utf-8")

def cat_nofilter_samples_percent(exp,sep,out):
    sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[sep]
    table = pd.DataFrame(pd.read_table(exp, sep=sep, index_col=0))
    table.columns = [str(i) for i in table.columns]
    table_name = table.index.name
    table = table.loc[(table > 0).any(axis=1), :]  # 去除全为0的行
    abund = table
    abund['Col_sum'] = abund.apply(lambda x: x.sum(), axis=1)
    abund_table = abund.sort_values(by=['Col_sum'], ascending=0)
    del abund_table["Col_sum"]  # 这两部是根据总量进行排序
    abund_table_path = os.path.join("taxa.table.xls")
    abund_table.to_csv(abund_table_path, sep="\t", encoding="utf-8")
    abund_table.columns = [str(i) for i in abund_table.columns]
    abund_table.loc['Row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
    sample_empt = []
    b = abund_table.apply(lambda x: x.sum(), axis=0)
    for i in range(len(abund_table.columns)):
        if b[i] > 0:
            pass
        else:
            sample_name = abund_table.columns[i]
            sample_empt.append(sample_name)
    abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['Row_sum'], axis=1).drop('Row_sum')
    abund_table_percent.to_csv(out, sep="\t", encoding="utf-8")


def get_others(per_exp,filter_value,out):
    perdf = pd.DataFrame(pd.read_table(per_exp, sep='\t', index_col=0))
    new_df = perdf.loc[list((perdf >= filter_value).any(axis=1)),]
    new_df2 = new_df.copy()
    others = perdf.loc[list((perdf < filter_value).all(axis=1))]
    if len(others) > 0:
        new_df2.loc["others"] = others.apply(lambda x: x.sum(), axis=0)
    other = os.path.join( "taxa.percents.table.xls")
    new_df2.to_csv(out, sep="\t", encoding="utf-8")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='\n该脚本用于根据传入的数据表格来获取数据绘制碱基柱状堆积图，需额外输入数据表格',
                                    add_help = False, usage = '\npython stacked.py -i exp_matrix -v value -o outfile')
    parser.add_argument('-i', type=str, metavar="input_matrix", required=True,help="please input exp file ")
    parser.add_argument('-s', type=str, metavar="sep_type", default="tab",required=True,help="please input sep type ")
    parser.add_argument('-v', type=float, metavar="threshole value", default=0.05, help="value for identiy others ", required=False)
    parser.add_argument('-o', type=str, metavar="output_dir",default=None, help="default is local dir. Output directory name", required=True)
    #
    args = parser.parse_args()

    input_exp,filter,out_result ,sep= args.i,args.v,args.o,args.s
    if filter == 0.0 or  filter == "" or  filter == 0 :
        cat_nofilter_samples_percent(input_exp,sep,out_result)
    else:
        cat_samples_percent(input_exp,sep)
        get_others("taxa.percents.table.xls",filter,out_result)

