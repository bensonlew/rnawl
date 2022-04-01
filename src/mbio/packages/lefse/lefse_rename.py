# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import re,argparse

def file_biom(lefse,namefile,output):
    file = open(output,'w')
    lefse_dict={}
    with open(lefse,'r') as f:
         lines = f.readlines()
         for line in lines:
             line = line.strip().split('\t')
             line2='\t'.join(line[1:])
             lefse_dict[line[0]]=line2
    rename_dict={}
    with open(namefile,'r') as f:
         lines = f.readlines()
         for line in lines:
             line = line.strip().split('\t')
             rename_dict[line[1]]=line[0]
    for name in set(lefse_dict).intersection(rename_dict):
        line =rename_dict[name]+'\t'+lefse_dict[name]
        file.write(line+'\n')


if __name__=='__main__':
    pars=argparse.ArgumentParser()
    pars.add_argument('-i',metavar='[lefse_lda.xls]')
    pars.add_argument('-s',metavar='[newname_file.xls]')
    pars.add_argument('-o',metavar='[OUT FILE]')
    args=pars.parse_args()
    file_biom(args.i,args.s,args.o)
