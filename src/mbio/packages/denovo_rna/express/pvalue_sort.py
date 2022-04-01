# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

"""无参转录组go/kegg富集调控分析排序"""

import os
import re
import types


def go_sort_pval(enrich_path, regulate_path, stat_path, pvalue_path):
    with open(enrich_path, 'rb') as f1, open(regulate_path, 'rb') as f2, open(stat_path, 'wb') as w:
        w.write('Term_type\tTerm\tgo_id\tup_num\tup_percent\tdown_num\tdown_percent\tcorrected_pvalue\tpvalue\n')
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        for line1 in lines1[1:]:
            line1 = line1.strip().split('\t')
            for line2 in lines2:
                line2 = line2.strip().split('\t')
                if line1[0] == line2[2]:
                    line = line2[:7] + line1[6] + line1[-2]
                    for i in range(len(line)):
                        w.write(line[i] + '\t')
                    w.write('\n')
    with open(stat_path, 'rb') as f, open(pvalue_path, 'wb') as w:
        lines = f.readlines()
        for i in range(1, len(lines)):
            for j in range(1, len(lines)-1):
                m = lines[j].strip().split('\t')[-1]
                n = lines[j+1].strip().split('\t')[-1]
                if float(m) > float(n):
                    lines[j], lines[j+1] = lines[j+1], lines[j]
        for line in lines:
            w.write(line)
    return pvalue_path

def kegg_sort_pval(enrich_path, regulate_path, stat_path, pvalue_path):
    with open(enrich_path, "rb") as f, open(regulate_path, "rb") as r, open(stat_path, "wb") as w:
        w.write("term\tdatabase\tpathway id\tpvalue\tcorrected pvalue\tup_num\tdown_num\tup_percent\tdown_percent\n")
        lines1 = f.readlines()
        lines2 = r.readlines()
        for line1 in lines1:
            if re.match(r"\w", line1):
                line1 = line1.strip().split('\t')
                for line2 in lines2:
                    line2 = line2.strip().split('\t')
                    if line1[2] == line2[0]:
                        up_percent = float(float(line2[2])/float(line1[3]))
                        down_percent = float(float(line2[3])/float(line1[3]))
                        w.write(line1[0] + '\t' + line1[1] + '\t' + line1[2] + '\t' + line1[5] + '\t' + line1[6] + '\t' + line2[2] + '\t' + line2[3] + '\t')
                        w.write(str(up_percent) + '\t' + str(down_percent) + '\n')
    with open(stat_path, 'rb') as f, open(pvalue_path, 'wb') as w:
        lines = f.readlines()
        for i in range(1, len(lines)):
            for j in range(1, len(lines)-1):
                m = lines[j].strip().split('\t')[3]
                n = lines[j+1].strip().split('\t')[3]
                if float(m) > float(n):
                    lines[j], lines[j+1] = lines[j+1], lines[j]
        for line in lines:
            w.write(line)
    return pvalue_path
