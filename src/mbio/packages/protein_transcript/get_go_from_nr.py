#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import time
import os
from collections import defaultdict
import argparse


class GetGoFromNr(Controller):
    def __init__(self, idmapping, fasta_vs_nr, go_identity, out):
        super(GetGoFromNr, self).__init__()
        self.idmapping = os.path.abspath(idmapping)
        if not os.path.exists(self.idmapping):
            super(GetGoFromNr, self).end(normal=1, out='你传入的idmapping文件路径不正确')
        self.fasta_vs_nr = os.path.abspath(fasta_vs_nr)
        if not os.path.exists(self.fasta_vs_nr):
            super(GetGoFromNr, self).end(normal=1, out='你传入的nr比对结果文件路径不正确')
        self.out = os.path.abspath(out)
        if not os.path.exists(os.path.dirname(self.out)):
            super(GetGoFromNr, self).end(normal=1, out='输出文件夹路径不对')
        self.mapping_list = list()
        self.go_identity = go_identity

    def add_commands(self):
        with open(self.idmapping, 'r') as i_r:
            c = 0
            for n, line in enumerate(i_r):
                if n % 2000000 == 0:
                    c += 1
                    if not n == 0:
                        map_w.close()
                    mapping_c = self.add_command(name = 'idmapping' + str(c))
                    self.mapping_list.append(mapping_c)
                    map_w = open(os.path.join(mapping_c.work_dir, 'idmaping_splited.tab'), 'w')
                map_w.write(line)

    def how_run(self):
        self.after_some(self.mapping_list, self.set_output)

    def end(self):
        super(GetGoFromNr, self).end(out='your program finished')

    def run_mapping(self):
        for map in self.mapping_list:
            cmd = "python /mnt/ilustre/users/yitong.feng/scripts/annotion/get_GO_from_blast_by_nr.py %s %s %s %s"%(self.fasta_vs_nr, os.path.join(map.work_dir, 'idmaping_splited.tab'), self.go_identity, os.path.join(map.work_dir, 'go_list'))
            cmd += '\n' + 'rm ' + os.path.join(map.work_dir, 'idmaping_splited.tab') + '\n'
            params = dict(
                cmd=cmd,
                node=10,
                memory=50
            )
            map.set_params(params)
            map.run()

    def set_output(self):
        acc2go = defaultdict(list)
        for map in self.mapping_list:
            go_list = os.path.join(map.work_dir, 'go_list')
            if not os.path.exists(go_list):
                super(GetGoFromNr, self).end(normal=1, out='%s没有正确生成go list文件'%map.name)
            with open(go_list, 'r') as gr:
                for line in gr:
                    line = line.strip().split('\t')
                    if len(line) > 1:
                        acc2go[line[0]] += line[1].split('; ')
            with open(self.out, 'w') as gw:
                for acc in acc2go:
                    gw.write(acc + '\t' + '; '.join(acc2go[acc]) + '\n')
        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_mapping()
        self.fire()

if __name__ == '__main__':
    # def __init__(self, idmapping, fasta_vs_nr, go_identity, out):
    parser = argparse.ArgumentParser(description="The script to get go from nr")
    parser.add_argument("-idmapping", type=str, required=True, help="the idmapping tab")
    parser.add_argument("-fasta_vs_nr", type=str, required=True, help="the nr blast result file")
    parser.add_argument("-go_identity", type=str, default='98',help='for filter')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'nr.go.list'))

    args = parser.parse_args()
    mapping = GetGoFromNr(args.idmapping, args.fasta_vs_nr, args.go_identity, args.out)
    mapping.run()