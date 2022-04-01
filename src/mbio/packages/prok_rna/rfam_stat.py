# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
import re
import argparse


parser = argparse.ArgumentParser(description='small rna rfam数据库比对统计')
parser.add_argument('-i', '--xls', help='input rfam blast table', required=True)
parser.add_argument('-o', '--detail', help='output rfam blast table', required=True)
parser.add_argument('-db', '--rfam_seed', help='database of Rfam.seed', required=True)

args = vars(parser.parse_args())
blast = args['xls']
detail = args['detail']
rfam_seed = args['rfam_seed']

seed = {}
rfam_ac = ""
rfam_id = ""
rfam_de = ""
rfam_tp = ""
with open(rfam_seed, "r") as f:
    lines = f.readlines()
    for line in lines:
        m = re.match(r"#=GF\s+AC\s+(\S+)\s*$", line)
        m1 = re.match(r"#=GF\s+ID\s+(.*)\s*$", line)
        m2 = re.match(r"#=GF\s+DE\s+(.*)\s*$", line)
        m3 = re.match(r"#=GF\s+TP\s+(.*)\s*$", line)
        if m:
            rfam_ac = m.group(1)
        if m1:
            rfam_id = m1.group(1)
        if m2:
            rfam_de = m2.group(1)
        if m3:
            rfam_tp = m3.group(1)
            seed.update({rfam_ac:{"ID": rfam_id, "DE": rfam_de, "TP": rfam_tp}})
query = {}
columns = []
with open(blast, "r") as f, open (detail ,"w") as w:
    head = f.readline()
    w.write("query_name\tHit\tEvalue\tIdentity\tAC\tID\tTP\tDE\n")
    for line in f:
        items = line.strip().split("\t")
        rfam_ac = items[10].split(";")[0]
        hit = items[10].split(";")[1]
        if (not query.has_key(items[5])) and seed.has_key(rfam_ac):
            query[items[5]] = 1
            tp = seed[rfam_ac]["TP"]
            if re.match(r'\S+\s+rRNA.*', tp):
                columns = [items[5], hit, items[1], items[3], rfam_ac, seed[rfam_ac]["ID"], "rRNA", seed[rfam_ac]["DE"]]
                w.write("\t".join(columns) + "\n")
            if re.match(r'\S+\s+tRNA.*', tp):
                columns = [items[5], hit, items[1], items[3], rfam_ac, seed[rfam_ac]["ID"], "tRNA", seed[rfam_ac]["DE"]]
                w.write("\t".join(columns) + "\n")
            if re.match(r'\S+\s+snoRNA.*', tp):
                columns = [items[5], hit, items[1], items[3], rfam_ac, seed[rfam_ac]["ID"], "snoRNA", seed[rfam_ac]["DE"]]
                w.write("\t".join(columns) + "\n")
            if re.match(r'\S+\s+snRNA.*', tp):
                columns = [items[5], hit, items[1], items[3], rfam_ac, seed[rfam_ac]["ID"], "snRNA", seed[rfam_ac]["DE"]]
                w.write("\t".join(columns) + "\n")
            # if re.match(r'\S+\s+sRNA.*', tp):
            #     columns = [items[5], hit, items[1], items[3], rfam_ac, seed[rfam_ac]["ID"], "other_sRNA", seed[rfam_ac]["DE"]]
            #     w.write("\t".join(columns) + "\n")
            else:
                pass

