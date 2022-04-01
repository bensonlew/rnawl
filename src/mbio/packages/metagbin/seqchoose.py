# -*- coding: utf-8 -*-
# author = xieshichang
# version = 0.1
import re
import os
import sys
import argparse
import time


def set_index(infile, idxfile, tag='>'):
    pre_seek = 0
    pre_name = 'seq_header'
    l_tag = len(tag)
    w = open(idxfile, 'w')
    r = open(infile, 'r')
    l = 'fake l'
    while l:
        l = r.readline()
        if l.startswith(tag):
            this_seek = r.tell() - len(l)
            this_name = l[l_tag:].strip()
            pre_len = this_seek - pre_seek
            w.write('{}\t{}\t{}\n'.format(pre_name, pre_seek, pre_len))
            pre_seek = this_seek
            pre_name = this_name
    w.write('{}\t{}\t{}\n'.format(pre_name, pre_seek, -1))
    r.close()
    w.close()


def extract(infile, idxfile, listfile, outfile):
    li = {}
    with open(listfile, 'rb') as lf:
        li = {l.strip():1 for l in lf}
    with open(infile, 'rb') as fin, open(idxfile, 'rb') as fidx, open(outfile, 'w') as w:
        for l in fidx:
            line = l.strip().split()
            if line[0] in li:
                fin.seek(int(line[-2]))
                w.write(fin.read(int(line[-1])))


def _main(args):
    if args.fa:
        tag = '>'
    if args.fq:
        tag = '@'
    tag = args.tag or tag

    idxfile = args.infile + '.idx'
    if os.path.exists(idxfile):
        os.remove(idxfile)
    t = time.time()
    print('seq choosing start')
    if not os.path.exists(idxfile):
        print('start indexing')
        set_index(args.infile, idxfile, tag)
        print('\tindexing done! time used: {}s'.format(time.time() - t))
    if args.outfile:
        t2 = time.time()
        print('start extract')
        extract(args.infile, idxfile, args.listfile, args.outfile)
        print('\textract done! time used: {}s'.format(time.time() - t2))
    print('seq choosing done! time used: {}s'.format(time.time() - t))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--infile', required=True)
    parser.add_argument('-l', '--listfile')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-t', '--tag', default=None)

    seq_type = parser.add_mutually_exclusive_group()
    seq_type.add_argument('-fa', action='store_true', default=True)
    seq_type.add_argument('-fq', action='store_true', default=False)

    args = parser.parse_args()
    _main(args)

