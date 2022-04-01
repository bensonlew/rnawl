#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__ = 'Wanjin.Hu'

import argparse
import os

parser = argparse.ArgumentParser(description="Rename reads in fastq files")
parser.add_argument("-i", "--input", dest="inFastq", required=True,
        type=str, help="input fastq file")
parser.add_argument("-o", "--output", dest="outFastq", required=True,
        type=str, help="output rename fastq file")
parser.add_argument("-n", "--name", dest="readsName", required=True,
        type=str, help="reads name which need rename")
args = parser.parse_args()

renameFq = open(args.outFastq,'w')
with open(args.inFastq,'r')as list:
    c = 1
    for i, line in enumerate(list):
        n = (i+4)/4
        wholeLine = line.strip('\n')
        if n % 1 == 0:
            str1 = wholeLine[0]
            str2 = str1+args.readsName+'_'+str(c)
            c = c + 1
        else:
            str2 = wholeLine
        renameFq.write(str2+'\n')
    renameFq.close()
