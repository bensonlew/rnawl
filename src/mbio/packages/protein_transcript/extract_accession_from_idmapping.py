#!/usr/bin/env python
#coding:utf-8
#fengyitong
#20190218

import sys
from multiprocessing import Pool


def get_acc_info(line):
    acc_list = list()
    if line.strip:
        items = line.strip().split('\t')
        acc = items[0]
        nr = items[3]
        if nr in nr_list:
            acc_list.append((acc,nr))
            nr_list.remove(nr)
    return acc_list

def chunk_file(in_file, chunk_size=100000):
    line_list = []
    i = 0
    with open(in_file, 'rb') as fr:
        for num, seq in enumerate(fr):
            line_list.append(seq)
            i += 1
            if i == chunk_size:
                yield line_list
                i = 0
                line_list = []
        else:
            if line_list:
                yield line_list

def deal_chunks(line_list):
    res = list()
    for line in line_list:
        res += get_acc_info(line)
    # map(get_conn_info, line_list)
    return res

def multi_process(in_file, p_num=8, chunk_size=100000):
    out_list = list()
    pool = Pool(p_num)
    results = list()
    for line_list in chunk_file(in_file, chunk_size):
        future = pool.apply_async(deal_chunks, (line_list,))
        while len(results) > 100:
            for fu in results:
                if fu.ready():
                    out_list += fu.get()
                    results.remove(fu)
        results.append(future)
    pool.close()
    pool.join()
    for fu in results:
        if fu.ready():
            out_list += fu.get()
    return out_list


if __name__ == '__main__':

    if len(sys.argv) != 4:
        exit('USAGE:\npython %s acc2nr.list idmappingDB out_file' % sys.argv[0])

    acc2nr = sys.argv[1]
    idmapping = sys.argv[2]
    out = sys.argv[3]

    with open(acc2nr, 'r') as ar:
        ar_info = ar.readlines()
        nr_list = [line.strip().split('\t')[1] for line in ar_info if line.strip()]
        nr2acc = {line.strip().split('\t')[1]:line.strip().split('\t')[0] for line in ar_info if line.strip()}

    out_list = multi_process(idmapping)

    with open(out, 'w') as out_w, open(out+'_relation', 'w') as out_c:
        # out_w.write('\n'.join(out_list))
        for out in out_list:
            out_w.write(out[0] + '\n')
            out_c.write(out[0] + '\t' + nr2acc[out[1]] + '\n')