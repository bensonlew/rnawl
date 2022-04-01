#!/usr/bin/env python
#coding:utf-8
#fengyitong
#20190218

import sys
# from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio.Blast import NCBIXML
import re

if len(sys.argv) != 6:
    exit('USAGE:\npython %s blast.list DB.txt delist identity out_file'%sys.argv[0])

# id2dbid=defaultdict(list)
id2dbid=dict()

blast=sys.argv[1]
db=sys.argv[2]
delist = sys.argv[3]
identity = float(sys.argv[4])
out=sys.argv[5]
out_detail='%s.detail.xls'%out

with open(delist, 'r') as de_r:
    des = [line.strip().split('\t')[0] for line in de_r if line.strip()]

if blast.endswith('.xls'):
    with open(blast, 'r') as blast_r:
        for line in blast_r.readlines():
            if line.strip():
                items = line.strip().split('\t')
                if items[0] in des and float(items[3]) >= identity:
                    if items[1] not in id2dbid:
                        id2dbid[items[1]] = list()
                    id2dbid[items[1]].append(items[2])
else:
    with open(blast, 'r') as blast_r:
        records = NCBIXML.parse(blast_r)
        for rec in records:
            query = re.split(' ', rec.query, maxsplit=1)[0]
            if query in des:
                for align in rec.alignments:
                    for hsp in align.hsps:
                        hit = align.hit_def
                        ident = float(hsp.identities)
                        hit_len = float(hsp.align_length)
                        if ident/hit_len*100 >= identity:
                            if hit not in id2dbid:
                                id2dbid[hit] = list()
                            id2dbid[hit].append(query)

accs = id2dbid.keys()
print(id2dbid)
out_list = list()

def get_conn_info(line):
    if line.strip:
        items = line.strip().split(' ')
        col1 = list()
        col2 = list()
        try:
            col1 = id2dbid[items[0]]
            col2 = id2dbid[items[1]]
        except:
            pass
        if col1 and col2:
            for i in col1:
                for j in col2:
                    out_list.append([i, j] + items)

db_r = open(db,'r')
header=db_r.readline().strip()
heads=header.strip().split(' ')
# with ThreadPoolExecutor(8) as pool:
#     pool.map(get_conn_info, db_r)
try:
    map(get_conn_info, db_r)
except:
    pass
db_r.close()

with open(out_detail,'w') as out_d, open(out, 'w') as out_w:
    out_w.write('protein1\tprotein2\tscore\n')
    out_d.write('protein1\tprotein2\t%s\n'%'\t'.join(heads))
    mark=0
    for i in range(len(heads)):
        if u'score' in heads[i]:
            mark=i
    for info in out_list:
        out_w.write('%s\t%s\t%s\n' % (info[0], info[1], info[mark + 2]))
        out_d.write('\t'.join(info) + '\n')