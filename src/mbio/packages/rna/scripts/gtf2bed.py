# -*- coding: utf-8 -*-
# __author__ = 'shijin'
import os
import re
import sys
import types
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="gtf转换bed文件")
parser.add_argument("-i", "--input", help="输入gtf文件")
parser.add_argument("-o", "--output", help="输出bed文件")

args = parser.parse_args()
if isinstance(args.input, types.NoneType):
    sys.exit(1)

chr = defaultdict(lambda: [])
dir = defaultdict(lambda: [])
exon_start = defaultdict(lambda: [])
exon_end = defaultdict(lambda: [])
cds_start = defaultdict(lambda: [])
cds_end = defaultdict(lambda: [])

r = open(args.input, "r")
for line in r:
    if line.startswith("#"):
        continue
    else:
        tmp = line.strip().split("\t")
        m = re.match("transcript_id \"(.+?)\"",tmp[-1])
        if m:
            transcript_id = m.group(1)
            chr[transcript_id].append(tmp[0])
            dir[transcript_id].append(tmp[6])
            if tmp[2] == "exon":
                exon_start[transcript_id].append(tmp[3])
                exon_end[transcript_id].append(tmp[4])
            elif tmp[2] == "CDS":
                cds_start[transcript_id].append(tmp[3])
                cds_end[transcript_id].append(tmp[4])
        else:
            print line
            break
                
f = open(args.output + "_tmp", "w")
for transcript_id in chr.keys():
    chr_lst = chr[transcript_id]
    dir_lst = dir[transcript_id]
    if len(set(chr_lst)) != 1 or len(set(dir_lst)) != 1:
        print "基因不在同一条染色体上或方向不一致，程序退出"
        print set(chr_lst)
        print set(dir)
        sys.exit(1)
    chr_num = chr_lst[0]  # 第一列
    start = str(min([int(x) for x in exon_start[transcript_id]]) - 1)  # 第二列
    end = str(max([int(x) for x in exon_end[transcript_id]]))  # 第三列, transcript_id为第四列
    five_col = "0"  # 第五列
    dirc = dir_lst[0]  # 第六列
    try:
        cds_st = str(min([int(x) for x in cds_start[transcript_id]]) - 1)  # 第七列
    except:
        cds_st = str(min([int(x) for x in exon_start[transcript_id]]) -1)
    try:
        cds_ed = str(max([int(x) for x in cds_end[transcript_id]]))  # 第八列，第九列为0
    except:
        cds_ed = str(max([int(x) for x in exon_end[transcript_id]]))
    ex_num = str(len(exon_start[transcript_id]))  # 第十列
    length_list = []
    for i in range(int(ex_num)):
        ex_start = int(exon_start[transcript_id][i])
        ex_end = int(exon_end[transcript_id][i])
        length = ex_end - ex_start + 1
        length_list.append(str(length))
    exsz = ",".join(length_list)  # 第十一列
    sta_list = []
    for i in range(int(ex_num)):
        if i == 0:
            mark = int(exon_start[transcript_id][i])
        ex_start = int(exon_start[transcript_id][i])
        sta_list.append(str(ex_start - mark))
    exst = ",".join(sta_list)  # 第十二列
    string = chr_num + "\t"
    string += start + "\t"
    string += end + "\t"
    string += transcript_id + "\t"
    string += five_col + "\t"
    string += dirc + "\t"
    string += cds_st + "\t"
    string += cds_ed + "\t"
    string += "0\t"
    string += ex_num + "\t"
    string += exsz + ",\t"
    string += exst + ",\n"
    f.write(string)
    
f.close()
r.close()
os.system("less {} | sort -n > {}".format(args.output + "_tmp", args.output))
    
    
    
