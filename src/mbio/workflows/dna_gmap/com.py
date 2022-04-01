#!/usr/bin/env python
# coding:utf-8
# fengyitong
# 20171122

import sys

if len(sys.argv) < 5:
    exit('python combined_2_id.py <idmap> <id1.list> <id2.list> <id1.new.list> <id2.new.list>\n')

from collections import defaultdict

idfile = sys.argv[1]
id1 = sys.argv[2]
id2 = sys.argv[3]
id3 = sys.argv[4]
id4 = sys.argv[5]

idmap1 = defaultdict(list)
idmap2 = defaultdict(list)

with open(idfile, 'r') as idfile_r:
    lines = idfile_r.readlines()
    for line in lines:
        ids = line.strip().split('\t')
        idmap1[ids[0]].append(ids[1])  # pro2rna
        idmap2[ids[1]].append(ids[0])  # rna2pro

with open(id1, 'r') as id1_r, open(id2, 'r') as id2_r, open(id3, 'w') as id3_w, open(id4, 'w') as id4_w:
    gene_list = [x.strip() for x in id1_r.readlines()]
    pro_list = [x.strip() for x in id2_r.readlines()]
    for pro in pro_list:
        if pro in idmap1:
            id3_w.write('%s|%s\n' % (pro, idmap1[pro][0]))
            id4_w.write('%s|%s\n' % (pro, idmap1[pro][0]))
            gene_list.remove(idmap1[pro][0])
        else:
            id3_w.write('%s\n' % (pro))
    for rna in gene_list:
        id4_w.write('%s\n' % (rna))
