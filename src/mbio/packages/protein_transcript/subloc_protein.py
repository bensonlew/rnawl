# -*- coding: utf-8 -*-
# author fengyitong 2019-02

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
import os
import argparse
from Bio import SeqIO
from collections import Counter


class SublocProtein(Controller):
    def __init__(self, fasta, specie, go_list, out):
        super(SublocProtein, self).__init__()
        os.system('export PATH="/mnt/ilustre/centos7users/yitong.feng/software/libsvm/libsvm-3.23:$PATH"')
        self.mulit = '/mnt/ilustre/centos7users/yitong.feng/software/MultiLoc2-26-10-2009/src/multiloc2_prediction.py'
        self.fasta = os.path.abspath(fasta)
        if not os.path.exists(self.fasta):
            super(SublocProtein, self).end(normal=1, out='fasta文件路径不正确')
        self.specie = specie
        self.go_list = os.path.abspath(go_list)
        if not os.path.exists(self.go_list):
            super(SublocProtein, self).end(normal=1, out='需要传入GO.list文件')
        self.out = os.path.abspath(out)
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(SublocProtein, self).end(normal=1, out='输出文件夹路径不对')
        self.sub_list = list()
        self.split_fasta = list()
        self.split_name = dict()
        self.run_split_fasta()

    def run_split_fasta(self):
        '''
        分割fasta序列
        '''
        line = 1
        i = 1
        w = open('fasta_%s' %i, 'wb')
        self.split_fasta.append('fasta_1')
        for seq_record in SeqIO.parse(self.fasta, "fasta"):
            if line <= 100:
                w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
                line += 1
            else:
                i += 1
                w.close()
                line = 1
                w = open('fasta_%s' % i, 'wb')
                w.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
                self.split_fasta.append('fasta_%s' %i)
                line += 1

            if self.split_name.has_key(i):
                self.split_name[i].append(seq_record.id)
            else:
                self.split_name[i] = [seq_record.id]
        w.close()

    def run_split_go(self):
        '''
        分割go注释文件
        '''
        go_dict = dict()
        with open(self.go_list, 'r') as go:
            for line in go.readlines():
                go_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
        for file_key in self.split_name.keys():
            genes = self.split_name[file_key]
            with open('go_%s' % file_key, 'wb') as go_split:
                for gene in genes:
                    if go_dict.has_key(gene):
                        go_split.write("{}\t{}\n".format(gene, go_dict[gene]))

    def add_commands(self):
        for i in self.split_name:
            sub = self.add_command(name='sub_%s'%i)
            self.sub_list.append(sub)

    def how_run(self):
        self.after_some(self.sub_list, self.set_output)

    def end(self):
        super(SublocProtein, self).end(out='your program finished')

    def run_subloc(self):
        self.run_split_go()
        db = {
            'Fungi': 'fungal',
            'Animals': 'animal',
            'Plants': 'plant',
        }
        if self.specie not in db:
            specie = 'fungal'
        else:
            specie = db[self.specie]
        for subloc in self.sub_list:
            file_key = subloc.name.split('sub_')[1]
            cmd = '{} {} '.format('python', self.mulit)
            cmd += '-{}={} '.format("origin", specie)
            cmd += '-{}={} '.format("fasta", '../fasta_%s' % file_key)
            cmd += '-{}={} '.format("go", '../go_%s' % file_key)
            cmd += '-{}={} '.format("result", 'multiloc_result')
            cmd += '-predictor=HighRes'
            params = dict(
                cmd=cmd,
                node=4,
                memory=12
            )
            subloc.set_params(params)
            subloc.run()

    def set_output(self):
        mul_rel = os.path.join(self.out, 'multiloc.xls')
        mul_stat = os.path.join(self.out, 'multiloc_stat.xls')
        first_loc = []
        for subloc in self.sub_list:
            with open(mul_rel, 'a') as multiloc_result:
                for file in os.listdir(subloc.work_dir):
                    if file.startswith('multiloc_result'):
                        with open(os.path.join(subloc.work_dir,file)) as multiloc:
                            lines = multiloc.readlines()
                            for line in lines:
                                if line.startswith(('MultiLoc2 ', 'predictor =', 'origin =')) or not line.strip():
                                    continue
                                line = line.strip().split('\t')
                                if len(line) > 8:
                                    multiloc_result.write('\t'.join(line[0:4]) + '\n')
                                    first_loc.append(line[1].split(':')[0])
        loc_stat = Counter(first_loc)
        stat_list = sorted(loc_stat.items(), key=lambda ot: ot[1], reverse=True)
        with open(mul_stat, 'w') as multiloc_stat:
            multiloc_stat.write('Subcelluar location\tProtein num' + '\n')
            for loc in stat_list:
                multiloc_stat.write(loc[0] + '\t' + str(loc[1]) + '\n')
        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_subloc()
        self.fire()

if __name__ == '__main__':
    # def __init__(self, fasta, denovo, species, go_list, out)
    parser = argparse.ArgumentParser(description="The script to run protein subloc analyse")
    parser.add_argument("-fasta", type=str, required=True, help="the protein fasta")
    parser.add_argument("-specie", type=str, default='Animals', help="Animals,Plants,Fungi,you can choose one")
    parser.add_argument("-go_list", type=str, default='no',help='your GO.list')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'subloc_results'))

    args = parser.parse_args()
    subloc = SublocProtein(args.fasta, args.specie, args.go_list, args.out)
    subloc.run()