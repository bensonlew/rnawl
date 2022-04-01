#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# __version__ = 'v1.0'
# __last_modified__ = '20151120'
"""
将七个水平的otu_taxon_table整合成lefse的输入文件的格式
"""

import os
import argparse

def get_argu():
    par = argparse.ArgumentParser()
    par.add_argument("-i", metavar="[taxon_sum_dir]", required=True, help="输入tax_summary_a文件夹")
    par.add_argument("-g", metavar="[groupfile]", required=True, help="输入分组文件信息")
    par.add_argument("-o", metavar="[outfile]", required=True, help="输入输出文件名")
    args = par.parse_args()
    return args

def get_index(choose, oldlist):
    newindex = [oldlist.index(i) for i in choose]
    return newindex

def get_lefse_input(gfile, output, taxondir):
    a = open(gfile)
    lefse_input = open(output, 'w')
    alines = a.readlines()
    groupnum = len(alines[0].split()) - 1
    sample = []
    group1 = []
    if groupnum == 1:
        for i in alines:
            sample.append(i.split()[0])
            group1.append(i.split()[1])
        lefse_input.write('\t'.join(group1) + '\n')
        lefse_input.write('\t'.join(sample) + '\n')
    else:
        group2 = []
        for i in alines:
            sample.append(i.split()[0])
            group1.append(i.split()[1])
            group2.append(i.split()[2])
        lefse_input.write('\t'.join(group1) + '\n')
        lefse_input.write('\t'.join(group2) + '\n')
        lefse_input.write('\t'.join(sample) + '\n')
    def get_values(myfile):
        temp = open(myfile)
        templines = temp.readlines()
        index = get_index(sample[1:], templines[0].split()[1:])
        allnewline = []
        for i in templines[1:]:
            values = i.split()
            newtax = '|'.join(values[0].split(';'))
            myvalue = [values[1:][n-1] for n in index]
            myvalue.insert(0, newtax)
            allnewline.append(myvalue)
        temp.close()
        return allnewline

    files = []
    for i in os.listdir(taxondir):
        if '.txt' in i:
            files.append(taxondir+'/' + i)
    for i in files:
        allline = get_values(i)
        for line in allline:
            lefse_input.write('\t'.join(line) + '\n')
    a.close()
    lefse_input.close()

    #if groupnum == 1:
    #    os.system('/mnt/ilustre/users/sanger/app/meta/lefse/format_input.py  lefse_input.txt  lefse_format.txt  -f  r -c 1 -u 2 -o 1000000')
    #else:
    #    os.system('/mnt/ilustre/users/sanger/app/meta/lefse/format_input.py  lefse_input.txt  lefse_format.txt  -f  r -c 1 -s 2 -u 3 -o 1000000')


def main():
    opts = get_argu()
    get_lefse_input(opts.g, opts.o, opts.i)
if __name__ == '__main__':
    main()
