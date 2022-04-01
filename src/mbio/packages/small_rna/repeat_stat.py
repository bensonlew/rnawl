# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
import re
import argparse


parser = argparse.ArgumentParser(description='repeat比对统计')
parser.add_argument('-i', '--xls', help='input repeat blast table', required=True)
parser.add_argument('-o', '--detail', help='output repeat detail table', required=True)
parser.add_argument('-f', '--gff', help='repeat annotation file', required=True)

args = vars(parser.parse_args())
blast = args['xls']
detail = args['detail']
gff = args['gff']

ssr = {}
with open(gff, "r") as f:
    head = f.readline()
    for line in f:
        items = line.strip().split("\t")
        info = re.search(r'ID=(\S+);.*;Class=(\S+).*;PercDiv=.*', items[8])
        ID = info.group(1)
        CLASS = info.group(2).split("/")[0]
        ssr[ID] = CLASS

query = {}
columns = []
with open(blast, "r") as f, open (detail ,"w") as w:
    head = f.readline()
    w.write("query_name\tHit\tEvalue\tIdentity\tClass\n")
    for line in f:
        items = line.strip().split("\t")
        if not query.has_key(items[5]):
            query[items[5]] = 1
            columns = [items[5], items[10], items[1], items[3], ssr[items[10]]]
            w.write("\t".join(columns) + "\n")
        else:
            pass