# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
import os
import re
import argparse
from collections import defaultdict


class RunAnnotation(Controller):
    def __init__(self, exp_fa, kegg_db, go_db, string_db, kegg_evalue, go_evalue, cog_evalue, out):
        super(RunAnnotation, self).__init__()
        self.exp_fa = exp_fa
        if not os.path.exists(self.exp_fa):
            # exit('你传入的fasta文件不存在')
            super(RunAnnotation, self).end(normal=1,out='你传入的fasta文件不存在')
        self.kegg_db = kegg_db
        if not os.path.exists(self.kegg_db):
            super(RunAnnotation, self).end(normal=1, out='你传入的keggdb文件不存在')
        self.build_kegg = 0
        if not self.kegg_db.endswith('.dmnd'):
            self.build_kegg = 1
        self.go_db = go_db
        if not os.path.exists(self.go_db):
            super(RunAnnotation, self).end(normal=1, out='你传入的godb文件不存在')
        self.string_db = string_db
        if not os.path.exists(self.string_db):
            super(RunAnnotation, self).end(normal=1, out='你传入的stringdb文件不存在')
        self.kegg_evalue = kegg_evalue
        self.go_evalue = go_evalue
        self.cog_evalue = cog_evalue
        self.out = out
        self.fsp_list = self.split_fasta(self.exp_fa)
        self.kegg_commands = list()
        self.go_commands = list()
        self.string_commands = list()
        self.can_end = 0

    def split_fasta(self, fa):
        sp_fas = list()
        n = 1
        na2fa = defaultdict(int)
        with open(fa, 'r') as fa_r:
            fa_info = fa_r.read().lstrip('>').split('\n>')
            limit_num = len(fa_info)/10
            for i in range(len(fa_info)):
                with open('splited' + str(n) + '.fa', 'a') as sp_w:
                    if n not in na2fa:
                        na2fa[n] = 0
                    if na2fa[n] == int(limit_num) - 1:
                        sp_w.write('>' + fa_info[i] + '\n')
                        na2fa[n] += 1
                        sp_fas.append('splited' + str(n) + '.fa')
                        n += 1
                    else:
                        sp_w.write('>' + fa_info[i] + '\n')
                        na2fa[n] += 1
            if 'splited' + str(n) + '.fa' not in sp_fas and os.path.exists('splited' + str(n) + '.fa'):
                sp_fas.append('splited' + str(n) + '.fa')
            return sp_fas

    def add_commands(self):
        if self.build_kegg:
            self.build_kegg_c = self.add_command(name='build_kegg')
        for fa in self.fsp_list:
            c_go = self.add_command(name=fa+'_go')
            c_kegg = self.add_command(name=fa+'_kegg')
            c_string = self.add_command(name=fa+'_string')
            self.go_commands.append(c_go)
            self.kegg_commands.append(c_kegg)
            self.string_commands.append(c_string)
        self.string2cog = self.add_command(name='string2cog')

    def how_run(self):
        if self.build_kegg:
            self.after_one(self.build_kegg, self.run_kegg)
        self.after_some(self.go_commands, self.merge_go)
        self.after_some(self.kegg_commands, self.merge_kegg)
        self.after_some(self.string_commands, self.merge_string)
        end_list = self.go_commands + self.kegg_commands
        end_list.append(self.string2cog)
        self.after_some(end_list, self.end)

    def end(self):
        print(self.can_end)
        if self.can_end == 3:
            super(RunAnnotation, self).end(out='your program finished')
        else:
            self.init_func(self.end)

    def build_kegg_db(self):
        db = os.path.join(self.build_kegg_c.work_dir, os.path.basename(self.kegg_db))
        cmd = 'diamond makedb --in %s --db %s --threads 10'%(self.kegg_db, db)
        params = dict(
            cmd=cmd,
            node=10,
            memory=50
        )
        self.build_kegg_c.set_params(params)
        self.build_kegg_c.run()
        self.kegg_db = db

    def run_kegg(self):
        for c in self.kegg_commands:
            name = c.name.split('_')[0]
            fa = os.path.abspath(name)
            cmd = "diamond blastp -q %s -d %s --more-sensitive -e %s -o %s_vs_KEGG.xls -f 6 --threads 10" %(fa, self.kegg_db,self.kegg_evalue, name)
            params = dict(
                cmd=cmd,
                node=10,
                memory=50
            )
            c.set_params(params)
            c.run()

    def run_go(self):
        for c in self.go_commands:
            name = c.name.split('_')[0]
            fa = os.path.abspath(name)
            cmd = "diamond blastp -q %s -d %s --more-sensitive -e %s -o %s_vs_nr.fast.xls -f 6 --threads 10" %(fa, self.go_db,self.go_evalue, name)
            params = dict(
                cmd=cmd,
                node=10,
                memory=50
            )
            c.set_params(params)
            c.run()

    def run_string(self):
        for c in self.string_commands:
            name = c.name.split('_')[0]
            fa = os.path.abspath(name)
            cmd = "diamond blastp -q %s -d %s --more-sensitive -e %s -o %s_vs_string.xml -f 5 --threads 10" %(fa, self.string_db,self.cog_evalue, name)
            params = dict(
                cmd=cmd,
                node=10,
                memory=50
            )
            c.set_params(params)
            c.run()

    def merge_go(self):
        file_list = list()
        for c in self.go_commands:
            name = c.name.split('_')[0]
            file = os.path.join(c.work_dir, name + '_vs_nr.fast.xls')
            if not os.path.exists(file):
                # exit('the output of %s does not exists'%c.name)
                super(RunAnnotation, self).end(normal=1, out='the output of %s does not exists'%c.name)
            file_list.append(file)
        self.merge_diomand_out(file_list, type_='xls', out = os.path.join(self.out, 'exp.fasta_vs_nr.fast.xls'))
        self.can_end += 1
        print('go_end')

    def merge_kegg(self):
        file_list = list()
        out = os.path.join(self.out, 'exp.fasta_vs_KEGG.xls')
        for c in self.kegg_commands:
            name = c.name.split('_')[0]
            file = os.path.join(c.work_dir, name + '_vs_KEGG.xls')
            if not os.path.exists(file):
                super(RunAnnotation, self).end(normal=1, out='the output of %s does not exists' % c.name)
            file_list.append(file)
        self.merge_diomand_out(file_list, type_='xls', out = out)
        if self.kegg_db.endswith('ko.pep.dmnd'):
            self.change_kegg(out)
        else:
            self.can_end += 1
            print('kegg_end')

    def merge_string(self):
        file_list = list()
        for c in self.string_commands:
            name = c.name.split('_')[0]
            file = os.path.join(c.work_dir, name + '_vs_string.xml')
            if not os.path.exists(file):
                super(RunAnnotation, self).end(normal=1, out='the output of %s does not exists' % c.name)
            file_list.append(file)
        self.merge_diomand_out(file_list, type_='xml', out = os.path.join(self.out, 'exp.fasta_vs_string.xml'))
        self.run_string2cog()
        self.can_end += 1
        print('string_end')

    def run_string2cog(self):
        cmd = 'perl /mnt/ilustre/users/bingxu.liu/workspace/annotation/String2Cog.pl -i {} --format blastxml -db {} -e 1e-3 -o {}/tmp_out'.format(os.path.abspath('exp.fasta_vs_string.xml'), '/mnt/ilustre/users/yitong.feng/scripts/annotion/cog.db', self.out)
        params = dict(
            cmd=cmd,
            node=10,
            memory=50
        )
        self.string2cog.set_params(params)
        self.string2cog.run()

    def change_kegg(self, kegg_out):
        with open('/mnt/ilustre/users/yitong.feng/scripts/annotion/ko_genes.list', 'r') as kg:
            kg2ko = {line.strip().split('\t')[1]:line.strip().split('\t')[0].split('ko:')[1] for line in kg if line.strip()}
        with open(kegg_out, 'r') as kor:
            kegg_info = kor.readlines()
        with open(kegg_out, 'w') as kow:
            for line in kegg_info:
                line = line.strip().split('\t')
                try:
                    tmp = line[1]
                    line[1] = kg2ko[tmp]
                    kow.write('\t'.join(line) + '\n')
                except:
                    pass
        self.can_end += 1
        print('kegg_end')

    def merge_diomand_out(self, file_list, type_ = 'xml', out = None):
        if type_ == 'xml':
            w = open(out, 'wb')
            head = False
            for f in file_list:
                with open(f, 'rb') as r:
                    flag = False
                    for line in r:
                        if not flag:
                            if not head:
                                w.write(line)
                            if re.match(r'^<Iteration>', line):
                                if head:
                                    w.write(line)
                                else:
                                    head = True
                                flag = True
                        elif not re.match(r'</BlastOutput', line) and line != '\n':
                            w.write(line)
                        else:
                            pass
            w.write(self.line_end)
            w.close()
        else:
            w = open(out, 'wb')
            head = True
            for f in file_list:
                with open(f, 'rb') as r:
                    for line in r:
                        if not head:
                            w.write(line)
                            head = True
                        else:
                            w.write(line)
            w.close()

    def run(self):
        self.add_commands()
        self.how_run()
        if self.build_kegg:
            self.build_kegg_db()
        else:
            self.run_kegg()
        self.run_go()
        self.run_string()
        super(RunAnnotation, self).fire()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="The pipeline for protein to run_annotation")
    parser.add_argument("-fasta", type=str, required=True, help="the query fasta")
    parser.add_argument("-go_db", type=str, required=True, help="go.dmnd")
    parser.add_argument("-kegg_db", type=str, required=True, help="kegg.dmnd or kegg.pep")
    parser.add_argument("-string_db", type=str, required=True, help="string.dmnd")
    parser.add_argument("-go_evalue", type=str, default='1e-5')
    parser.add_argument("-kegg_evalue", type=str, default='1e-5')
    parser.add_argument("-cog_evalue", type=str, default='1e-5')
    parser.add_argument("-out", type=str, default=os.getcwd())

    args = parser.parse_args()
    anno = RunAnnotation(args.fasta, args.kegg_db, args.go_db, args.string_db, args.kegg_evalue, args.go_evalue, args.cog_evalue, args.out)
    anno.run()