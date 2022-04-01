# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171214

"""统计比对到rfam数据库的xml文件的TP数"""
import os
import re
from Bio.Blast import NCBIXML
import argparse


parser = argparse.ArgumentParser(description='统计比对到rfam数据库的xml文件的TP数')
parser.add_argument('-i', '--xml', help='input xml file of blast rfam database', required=True)
parser.add_argument('-o', '--summary', help='output stat file', required=True)
parser.add_argument('-db', '--rfam_seed', help='database of Rfam.seed', required=True)

args = vars(parser.parse_args())
xml_fp = args['xml']
summary_out = args['summary']
rfam_seed = args['rfam_seed']

seed = {}
with open(rfam_seed, "r") as f:
    lines = f.readlines()
    for line in lines:
        m = re.match(r"#=GF\s+AC\s+(\S+)\s*$", line)
        m1 = re.match(r"#=GF\s+TP\s+(.*)\s*$", line)
        if m:
            rfam_id = m.group(1)
        if m1:
            seed[rfam_id] = m1.group(1)
rf_ids = {}
all_query = 0
unmatch = 0
with open(xml_fp) as f:
    records = NCBIXML.parse(f)
    for rec in records:
        query = re.split(' ', rec.query, maxsplit=1)[0]
        all_query += 1
        if rec.alignments:
            rf_id = rec.alignments[0].hit_id.split(";")[0]
            # rf_id = rec.alignments[0].hit_def.split(";")[0]
            if rf_id not in rf_ids.keys():
                rf_ids[rf_id] = 1
            else:
                rf_ids[rf_id] += 1
        else:
            unmatch += 1
taxon = {}
for rf_id in rf_ids.keys():
    if seed[rf_id] not in taxon.keys():
        taxon[seed[rf_id]] = rf_ids[rf_id]
    else:
        taxon[seed[rf_id]] += rf_ids[rf_id]
with open(summary_out, "w") as w:
    match = all_query - unmatch
    match_rate = round(float(match) / all_query * 100, 2)
    w.write("TP\ttotal_num\ttotal_percent\n")
    w.write("Total_all\t" + str(all_query) + "\t100%\n")
    w.write("All_matched\t" + str(match) + "\t" + str(match_rate) + "%\n")
    for ta in taxon.keys():
        rate = round(float(taxon[ta]) / all_query * 100, 2)
        w.write(ta + "\t" + str(taxon[ta]) + "\t" + str(rate) + "%\n")
