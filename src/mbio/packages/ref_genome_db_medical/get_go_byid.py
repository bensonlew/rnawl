#!/usr/bin/env python
#coding:utf-8
#liubinxu

import sys
import os
import datetime
from collections import defaultdict
from multiprocessing import Pool

acc_set = set()

def add_go(line):
    '''
    映射一条数据库NR与GO对应关系
    '''
    acc2go1 = []

    items = line.strip('\n').split('\t')
    acc_now = []
    if items[7] != "":
        for n,acc in enumerate(items):
            if n in [0,1,3,20,21] and acc != "":
                acc_now.extend([x.strip() for x in acc.split(";")])
        acc_now_set = set(acc_now)
        # print NR2GO().acc_set
        # exit
        global acc_set
        for acc in acc_now_set.intersection(acc_set):
            acc2go1.append("{}\t{}".format(acc, items[7]))
        return "\n".join(acc2go1)
    else:
        return ""

class ID2GO(object):
    def __init__(self):
        self.acc_set = ()
        self.idmapping_db = ""
        self.acc2go = []
        self.table = ""

    def get_acc_set(self, table):
        '''
        获取gene ID(NR比对结果)列表， swissprot 与数据库idmapping差异 去除 .1 .2
        '''
        # print "start acc"
        # print datetime.datetime.now().strftime('%H-%M-%S')
        self.table = table
        acc_list = []
        with open(table,'r') as fa_r:
            fa_r.readline()
            for line in fa_r:
                items = line.strip('\n').split('\t')
                acc_list.extend(items)
        self.acc_set = set(acc_list)
        global acc_set
        acc_set = self.acc_set
        # print "end acc"
        # print datetime.datetime.now().strftime('%H-%M-%S')

    def mapping_go(self, thread=10):
        '''
        映射nr_acc id  go对应关系
        '''
        # print "start mapping"
        # print datetime.datetime.now().strftime('%H-%M-%S')

        p = Pool(thread)
        idmapping_r = open(self.idmapping_db, 'r')
        self.acc2go = p.map(add_go, idmapping_r)
        p.close()
        p.join()
        idmapping_r.close()
        '''
        with open("nracc2go",'w') as go_w:
            go_w.write("*".join(self.acc2go))
        '''
        # print "end mapping"
        # print datetime.datetime.now().strftime('%H-%M-%S')

    def write_result(self, result_file):
        # print "start write result"
        # print datetime.datetime.now().strftime('%H-%M-%S')
        acc2go_dict = dict()

        for go_annot in self.acc2go:
            if go_annot == "":
                continue
            for nr2gos in go_annot.split("\n"):
                go_annot = nr2gos.split("\t")
                if acc2go_dict.has_key(go_annot[0]):
                    go = acc2go_dict[go_annot[0]]
                    go_new = [x.strip() for x in go_annot[1].split(";")]
                    acc2go_dict[go_annot[0]] = list(set(go + go_new))
                else:
                    acc2go_dict.update({go_annot[0]: [x.strip() for x in go_annot[1].split(";")]})

        with open(self.table,'r')as fa_r, open(result_file,'w')as nr_w:
            fa_r.readline()
            for line in fa_r:
                items = line.strip('\n').split('\t')
                gos = list()
                for item in items:
                    if acc2go_dict.has_key(item):
                        gos.extend(acc2go_dict[item])
                if len(gos) > 0:
                    nr_w.write("{}\t{}\n".format(items[0], ";".join(gos)))
        # print "end write result"
        # print datetime.datetime.now().strftime('%H-%M-%S')

if __name__ == "__main__":
    if len(sys.argv) - 1 != 4:
        exit('%s exp.fasta_vs_nr.xls idmapping.tb 10 nr.GO.list' %sys.argv[0])
    id2go = ID2GO()
    table = sys.argv[1]
    id2go.get_acc_set(table)
    id2go.idmapping_db = sys.argv[2]
    id2go.mapping_go(int(sys.argv[3]))
    id2go.write_result(sys.argv[4])
