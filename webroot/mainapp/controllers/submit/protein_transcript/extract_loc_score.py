# -*- coding: utf-8 -*-

import os
import sys
import re

if len(sys.argv) != 3:
    exit('USAGE: %s infile output' %sys.argv[0])


acc2loc = dict()
with open(sys.argv[1], 'r') as i_r:
    i_r.readline()
    for line in i_r:
        line = line.strip().split('\t')
        if len(line) >1 and ' ' in line[0]:
            acc_tmp = line[0].split(' ')[0]
            acc_list = list()
            start_list = list()
            for acc_loc in line[0].split('; '):
                if acc_loc.startswith('['):
                    acc_list.append(acc_tmp)
                    start_list.append(acc_loc.lstrip('[').split('-')[0])
                else:
                    acc_list.append(acc_loc.split(' ')[0])
                    start_list.append(acc_loc.split(' ')[1].lstrip('[').split('-')[0])
            acc2start = zip(acc_list, start_list)
            # acc = line[0].split(' ')[0]
            # accs = acc.split(';')
            # print(line[0].split(' ')[1])
            # start = line[0].split(' ')[1].lstrip('[').split('-')[0]
            for acc in acc2start:
                if acc[0] not in acc2loc:
                    acc2loc[acc[0]] = dict()
            for bl in line[1].split('; '):
                score = float(bl.split(': ')[1])
                if score < 75:
                    continue
                type_ = filter(str.isalpha, bl.split(': ')[0].split('(')[0])
                loc_ = filter(str.isdigit, bl.split(': ')[0].split('(')[0])
                for acc_ in acc2start:
                    acc, start = acc_
                    loc = int(start) + int(loc_) -1
                    if loc not in acc2loc[acc]:
                        acc2loc[acc][loc] = list()
                    tmp_dict = dict(type_=type_, score=score)
                    judge = 0
                    for i in acc2loc[acc][loc]:
                        if 'type_' in i and i['type_'] == type_:
                            judge = 1
                            if i['score'] > score:
                                continue
                            else:
                                i.update(**tmp_dict)
                    if not judge:
                        acc2loc[acc][loc].append(tmp_dict)
with open(sys.argv[2] + '1.txt', 'w') as w1, open(sys.argv[2] + '2.txt', 'w') as w2:
    for acc in acc2loc:
        for loc in acc2loc[acc]:
            tmp_list = [x['type_']+str(loc)+'('+str(x['score']) + ')' for x in acc2loc[acc][loc]]
            tmp_list2 = [acc + '\t' + x for x in tmp_list]
            w1.write(acc + '\t' + '\t'.join(tmp_list) + '\n')
            w2.write('\n'.join(tmp_list2) + '\n')
