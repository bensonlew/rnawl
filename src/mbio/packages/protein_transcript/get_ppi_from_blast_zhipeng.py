#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/21 16:55
"""
import sys
from Bio.Blast import NCBIXML
from queue import Queue
import multiprocessing
import time
import threading
import re

p_list = []


def split_file(file, com_queue, block_num):
    temp_list = []
    flag = 0
    with open(file) as infile:
        for line in infile:
            flag += 1
            temp_list.append(line)
            if flag == block_num:
                flag = 0
                com_queue.put(temp_list)
                temp_list = []
        if temp_list:
            com_queue.put(temp_list)
        com_queue.put(None)


def check_process(com_queue, process_num):
    while True:
        if len(p_list) < process_num:
            data = com_queue.get()
            if data is None:
                break
            p = multiprocessing.Process(target=p_run, args=(data,))
            p.start()
            p_list.append(p)
            continue
        while True:
            time.sleep(1)
            flag = False
            for p in p_list:
                if not p.is_alive():
                    p_list.remove(p)
                    flag = True
            if flag is True:
                break
    # 任务投递完成后join所有进程
    for p in p_list:
        p.join()

def p_run(data_list):
    '''处理代码'''
    print('===========================', multiprocessing.current_process().name)
    for line in data_list:
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
    print('+++++++++++++++++++++++++++', multiprocessing.current_process().name)


if __name__ == '__main__':
    if len(sys.argv) != 6:
        exit('USAGE:\npython %s blast.list DB.txt delist identity out_file' % sys.argv[0])

    # id2dbid=defaultdict(list)
    id2dbid = dict()

    blast = sys.argv[1]
    db = sys.argv[2]
    delist = sys.argv[3]
    identity = float(sys.argv[4])
    out = sys.argv[5]
    out_detail = '%s.detail.xls' % out

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
                            if ident / hit_len * 100 >= identity:
                                if hit not in id2dbid:
                                    id2dbid[hit] = list()
                                id2dbid[hit].append(query)
    db_r = open(db,'r')
    header=db_r.readline().strip()
    heads=header.strip().split(' ')

    s = time.time()
    # 并行数量
    parallel_num = 10
    # 块大小（行数，大文件建议 10万行）
    block_line_num = 100000
    queue = Queue(parallel_num + 5)
    out_list = list()
    thread1 = threading.Thread(target=split_file, args=(db, queue, block_line_num))
    thread2 = threading.Thread(target=check_process, args=(queue, parallel_num))
    thread1.start()
    thread2.start()
    thread1.join()
    thread2.join()
    print(time.time() - s)

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