#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/18 8:38
@file    : get_fasta_from_verylarge_db.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import sys
import gzip
from multiprocessing import Pool


def chunk_file(in_file, chunk_size=100000):
    id2seq = dict()
    if not in_file.endswith('.gz'):
        fr = open(in_file, 'r')
    else:
        fr = gzip.open(in_file, 'r')
    for line in fr:
        if len(id2seq) > chunk_size:
            yield id2seq
            id2seq = dict()
        if line.strip():
            if line.startswith('>'):
                # ids = [a.strip().lstrip('>').split(' ')[0] for a in line.strip().split('\x01')]
                ids = line.strip().lstrip('>').split('\x01')
                ids = [id.split(' ')[0] for id in ids]
                if u'|' in ids[0]:
                    ids = [id.split('|')[1] for id in ids]
                # if u'|' in id:
                #     id = id.split('|')[1]
            else:
                for id in ids:
                    tmp = ''
                    try:
                        tmp = id2seq[id]
                    except:
                        pass
                    if not tmp:
                        id2seq[id] = ''
                for id in ids:
                    id2seq[id] += line
    else:
        if id2seq:
            yield id2seq
    fr.close()

def filter_dict(raw_dict, exp_dict):
    filter_dict = dict()
    for exp in raw_dict:
        tmp = ''
        try:
            # tmp = exp_dict[exp]
            tmp = exp_dict[exp.split(' ')[0]]
        except:
            pass
        if not tmp:
            continue
        filter_dict[exp] = raw_dict[exp]
    # print(filter_dict)
    return filter_dict

def multi_process(in_file, exp_dict, p_num=8, chunk_size=100000):
    write_dict = dict()
    pool = Pool(p_num)
    results = list()
    for raw_dict in chunk_file(in_file, chunk_size):
        future = pool.apply_async(filter_dict, (raw_dict, exp_dict))
        while len(results) > 100:
            for fu in results:
                if fu.ready():
                    write_dict.update(fu.get())
                    results.remove(fu)
        results.append(future)
    pool.close()
    pool.join()
    for fu in results:
        if fu.ready():
            write_dict.update(fu.get())
    return write_dict


if __name__ == '__main__':

    if len(sys.argv) not in [4,5]:
        exit('USAGE:\npython %s exp.list DB.fasta chunk_size out_file' % sys.argv[0])

    exp = sys.argv[1]
    DB = sys.argv[2]
    if len(sys.argv) == 5:
        chunk_size = int(sys.argv[3])
        out_file = sys.argv[4]
    else:
        out_file = sys.argv[3]
        chunk_size = 100000

    with open(exp) as er:
        # exp_list = er.read().strip().split('\n')
        # while '' in exp_list:
        #     exp_list.remove('')
        exp_dict = {line.strip().split('\t')[0]: '1' for line in er if line.strip()}

    write_dict = multi_process(DB, exp_dict, p_num=8, chunk_size=chunk_size)
    with open(out_file, 'w') as fw:
        for id, seq in write_dict.items():
            fw.write('>' + id + '\n' + seq.strip() + '\n')

    with open(out_file + '_except.list', 'w') as ew:
        for i in exp_dict:
            tmp = ''
            try:
                tmp = write_dict[i]
            except:
                pass
            if not tmp:
                ew.write(i + '\n')