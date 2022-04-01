#!/bin/env python
# -*- coding: utf-8 -*- 
import argparse
import re
from mbio.packages.statistical.memory_profiler import profile
import pandas as pd


def get_argu():
    par = argparse.ArgumentParser()
    par.add_argument("-i", metavar="[level1234 file]", required=True, help="输入注释四个水平的文件")
    par.add_argument("-o", metavar="[outfile]", required=True, help="输出go注释表")
    args = par.parse_args()
    return args

@profile
def get_data(input,output):
    with open(input,'r') as p,open(output,'w') as g:
        lines=p.readlines()
        g.write(lines[0])
        for i in range(1,len(lines)):
            temp = lines[i].strip().split('\t')
            genes =temp[-1].split(";")
            for j in range(0,len(genes)):
                 g.write("\t".join(temp[:-1]) + "\t" + genes[j] + '\n')

@profile
def select_data(input,output):
    data=pd.read_table(input,sep='\t')    
    del data['Seq Number']
    del data['Percent']
    head = data.columns[:-1]
    d=dict(list(data.groupby(by=['Seq List'])))
    with open (output,'w') as f:
        f.write('#Query\t'+'\t'.join(head)+'\n')
        for k in d.keys():
            dict1 = {}
            for j in range(0,len(head)):
                dict1[head[j]] = ";" .join(set(list(d[k][head[j]])))
            list1 = []
            for p in head:
                list1.append(dict1[p])
            f.write(k + '\t' + '\t'.join(list1) + '\n')

if __name__ == '__main__':
    opts = get_argu()
    result = "result.xls"
    get_data(opts.i,result)
    select_data(result,opts.o)
