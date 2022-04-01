# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

#import pandas as pd
#import xml.etree.ElementTree as ET
import sys


def mergeGO(annoout, gofile):
    fi = open(annoout).read()
    allrecords = fi.split('\n')
    # create dictionary
    d = {}
    for item in allrecords:
        if item!='':
            l = item.split('\t')
            seq = l[0]
            seq=seq.split(' ')[0]
            GO_list = l[1].split(";")
            desc = l[2]
            if d.has_key(seq):
                for GO in GO_list:
                    if GO not in d[seq][0]:
                        d[seq][0].append(GO)
                d[seq][1].append(desc)
            else:
                golist = []
                for GO in GO_list:
                    if GO not in golist:
                        golist.append(GO)
                delist = []
                delist.append(desc)
                newlist = []
                newlist.append(golist)
                newlist.append(delist)
                d[seq] = newlist
    # create file without desc info
    keys = list(d)
    with open(gofile, 'w') as f:
        for dickey in keys:
            go_list = list(set(d[dickey][0]))
            go_list.sort()
            f.write(dickey + '\t' + ';'.join(go_list))
            if dickey != keys[-1]:
                f.write('\n')

mergeGO(sys.argv[1], sys.argv[2])
