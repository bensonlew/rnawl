#!/bin/env python
# -*- coding: utf-8 -*- 
import argparse
import re
from mbio.packages.statistical.memory_profiler import profile

def get_argu():
    par = argparse.ArgumentParser()
    par.add_argument("-a", metavar="[go_anno.xls]", required=True, help="输入go注释gene表")
    par.add_argument("-i", metavar="[level1234 file]", required=True, help="输入注释四个水平的文件")
    par.add_argument("-o", metavar="[outfile]", required=True, help="输出go注释表")
    args = par.parse_args()
    return args

@profile
def remove_gene(file1,file2,output):
    list1=[]
    with open (file1,'r') as f:
        lines =f.readlines()
        for line in lines[1:]:
            gene =line.strip('\n').split("\t")[0]
            list1.append(gene)
    with open (file2,'r') as p,open (output,'w') as w:
        lines =p.readlines()
        w.write(lines[0])
        for line in lines[1:]:
            temp =line.strip('\n').split("\t")
            arry =temp[-1].split(";")
            genes = []
            for ge in arry:
                genes.append(ge)
            temp2= list(set(list1) & set(genes))
            if len(temp2) >=1:
                des =";".join(temp2)
                temp[-1]=des
                w.write("\t".join(temp) + "\n")

if __name__ == '__main__':
    opts = get_argu()
    dict =remove_gene(opts.a,opts.i,opts.o)

