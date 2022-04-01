# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import re,argparse

def file_biom(otufile,output,output2):
    file = open(output,'w')
    file2 = open(output2, 'w')
    with open(otufile,'r') as f:
         lines = f.readlines()
         name = lines[0]
         names=name.strip().split('\t')
         names.append('taxonomy')
         new_name='\t'.join(names)
         num =0
         file.write(new_name+'\n')   
         for line in lines[1:]:
             num +=1
             name2 = 'tax'+str(num)
             line = line.strip().split('\t')
             liness =line[0]+'\t'+name2
             file2.write(liness+'\n')
             if num % 2 ==1:
                 line.append('d__aaa;k__ccc;p__'+name2)
                 line='\t'.join(line)
                 file.write(line+'\n')
             else:
                 line.append('d__bbb;k__ddd;p__'+name2)
                 line='\t'.join(line)
                 file.write(line+'\n')


if __name__=='__main__':
    pars=argparse.ArgumentParser()
    pars.add_argument('-i',metavar='[otu_table.xls]')
    pars.add_argument('-o',metavar='[OUT FILE]')
    pars.add_argument('-s',metavar='[OUT2 FILE2]')
    args=pars.parse_args()
    file_biom(args.i,args.o,args.s)
