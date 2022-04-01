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


class SublocProteinPipeline(Controller):
    def __init__(self, fasta, denovo, species, go_list, out, go_identity):
        super(SublocProteinPipeline, self).__init__()
        os.system('export PATH="/mnt/ilustre/centos7users/yitong.feng/software/libsvm/libsvm-3.23:$PATH"')
        self.parafly = '/mnt/ilustre/app/rna/assemble_mrna/trinityrnaseq_r20140413/trinity-plugins/parafly-r2013-01-21/bin/ParaFly'
        self.mulit = '/mnt/ilustre/centos7users/yitong.feng/software/MultiLoc2-26-10-2009/src/multiloc2_prediction.py'
        # self.mulit = '/mnt/ilustre/centos7users/yitong.feng/software/Multiloc_sanger/MultiLoc2-26-10-2009/src/multiloc2_prediction.py'
        self.fasta = os.path.abspath(fasta)
        if not os.path.exists(self.fasta):
            super(SublocProteinPipeline, self).end(normal=1, out='fasta文件路径不正确')
        self.denovo = denovo.lower()
        self.species = species.split(',')
        if self.denovo == 'no' and not os.path.exists(go_list):
            super(SublocProteinPipeline, self).end(normal=1, out='如果选择不从头开始做注释，则需要传入GO.list文件')
        self.go_list = ''
        if self.denovo == 'no':
            self.go_list = os.path.abspath(go_list)
        self.out = os.path.abspath(out)
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(SublocProteinPipeline, self).end(normal=1, out='输出文件夹路径不对')
        self.sub_list = list()
        self.split_fasta = list()
        self.split_name = dict()
        self.run_split_fasta()
        self.go_identity = go_identity

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
        if self.denovo != 'no':
            self.annotation = self.add_command(name='annotation')
            self.go_result_nr = self.add_command(name='go_result_nr')
            self.go_result = self.add_command(name='go_result')
        for specie in self.species:
            specie = self.add_command(name='sub_' + specie)
            self.sub_list.append(specie)

    def how_run(self):
        if self.denovo != 'no':
            self.after_one(self.annotation, self.run_go_result_nr)
            self.after_one(self.go_result_nr, self.run_go_result)
            self.after_one(self.go_result, self.run_subloc)
        self.after_some(self.sub_list, self.set_output)

    def end(self):
        super(SublocProteinPipeline, self).end(out='your program finished')

    def run_annotation(self):
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/annotion/run_annotation.py -fasta ' + self.fasta
        cmd += ' -out ' + os.path.join(self.annotation.work_dir, 'diomand_results')
        cmd += ' -run_what ' + 'go'
        params = dict(
            cmd=cmd,
            node=1,
            memory=5,
            qsub=0
        )
        self.annotation.set_params(params)
        self.annotation.run()

    def run_go_result_nr(self):
        vsnr = os.path.join(self.go_result_nr.work_dir, 'exp.fasta_vs_nr.fast.xls')
        try:
            os.link(os.path.join(self.annotation.work_dir, 'diomand_results', 'exp.fasta_vs_nr.fast.xls'), vsnr)
        except:
            if not os.path.exists(vsnr):
                super(SublocProteinPipeline, self).end(normal=1, out='注释没有正常完成')
        # cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/annotion/get_GO_from_blast_by_nr.py ' + vsnr + ' /mnt/ilustre/users//bing.yang/DB/GO/idmapping.tb ' + self.go_identity + ' nr.GO.list\n'
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/annotion/get_go_from_nr.py'
        cmd += ' -idmapping ' + '/mnt/ilustre/users//bing.yang/DB/GO/idmapping.tb'
        cmd += ' -fasta_vs_nr ' + 'exp.fasta_vs_nr.fast.xls'
        cmd += ' -go_identity ' + self.go_identity
        cmd += ' -out ' + os.path.join(self.go_result.work_dir, 'nr.GO.list') + '\n'
        params = dict(
            cmd=cmd,
            node=3,
            memory=9,
            qsub=0
        )
        self.go_result_nr.set_params(params)
        self.go_result_nr.run()

    def run_go_result(self):
        nrgo = os.path.join(self.go_result.work_dir, 'nr.GO.list')
        exp_list = os.path.join(self.go_result.work_dir, 'exp.list')
        with open(self.fasta, 'r') as fa, open(exp_list, 'w') as lw:
            acc_list = list()
            for line in fa:
                if line.startswith('>'):
                    acc_list.append(line.strip().lstrip('>'))
            lw.write('\n'.join(acc_list))
        if not os.path.exists(nrgo):
            super(SublocProteinPipeline, self).end(normal=1, out='GO报错，nr数据库mapping没有正常完成')
        if not os.path.exists(exp_list):
            super(SublocProteinPipeline, self).end(normal=1, out='没有exp.list文件')
        bin = "/mnt/ilustre/users/ting.kuang/LABELFREE/bin"
        uniprot2go = "/mnt/ilustre/users/ting.kuang/ITRAQ/db/GO/2018-11.GO.list"
        cmd = '{bin}/get_annot_list.pl exp.list {uniprot2go} pir.GO.list\n'.format(bin=bin,
                                                                                   uniprot2go=uniprot2go)
        cmd += '{bin}/merge_2tab_file.py nr.GO.list,pir.GO.list GO.list\n'.format(bin=bin)
        params = dict(
            cmd=cmd,
            node=4,
            memory=24
        )
        self.go_result.set_params(params)
        self.go_result.run()

    def run_subloc(self):
        if self.denovo != 'no':
            self.go_list = os.path.join(self.go_result.work_dir, 'GO.list')
            if not os.path.exists(self.go_list):
                super(SublocProteinPipeline, self).end(normal=1, out='GO.list没有正常生成')
        self.run_split_go()

        db = {
            'Fungi': 'fungal',
            'Animals': 'animal',
            'Plants': 'plant',
        }
        for subloc in self.sub_list:
            specie = subloc.name.split('sub_')[1]
            if specie not in db:
                specie = 'fungal'
            else:
                specie = db[specie]
            para_file = os.path.join(subloc.work_dir, 'parrel_cmd')
            with open(para_file, 'wb') as f:
                for file_key in self.split_name.keys():
                    cmd = '{} {} '.format('python', self.mulit)
                    cmd += '-{}={} '.format("origin", specie)
                    # cmd += '-{}={} '.format("output", "advanced")
                    cmd += '-{}={} '.format("fasta", '../fasta_%s' % file_key)
                    cmd += '-{}={} '.format("go", '../go_%s' % file_key)
                    cmd += '-{}={} '.format("result", 'multiloc_result_%s' % file_key)
                    cmd += '-predictor=HighRes'
                    f.write("{}\n".format(cmd))
            cmd = '{} '.format(self.parafly)
            cmd += '-{} {} '.format("c", "parrel_cmd")
            cmd += '-{} {} '.format("CPU", 20)
            cmd += '-v -shuffle'
            params = dict(
                cmd=cmd,
                node=20,
                memory=60
            )
            subloc.set_params(params)
            subloc.run()

    def set_output(self):
        for subloc in self.sub_list:
            specie = subloc.name.split('sub_')[1]
            mul_rel = os.path.join(self.out, '%s_multiloc.xls'%specie)
            mul_stat = os.path.join(self.out, '%s_multiloc_stat.xls'%specie)
            first_loc = []
            with open(mul_rel, 'w') as multiloc_result:
                for file in os.listdir(subloc.work_dir):
                    if file.startswith('multiloc_result_'):
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
        if self.denovo != 'no':
            self.run_annotation()
        else:
            self.run_subloc()
        self.fire()

if __name__ == '__main__':
    # def __init__(self, fasta, denovo, species, go_list, out)
    parser = argparse.ArgumentParser(description="The script to run protein subloc analyse")
    parser.add_argument("-fasta", type=str, required=True, help="the protein fasta")
    parser.add_argument("-denovo", type=str, default='yes', help="if you have GO.list choose no")
    parser.add_argument("-species", type=str, default='Animals', help="Animals,Plants,Fungi,you can choose all or single, split by ,")
    parser.add_argument("-go_list", type=str, default='no',help='if -denovo you choose no,please input your GO.list')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'subloc_results'))
    parser.add_argument("-go_identity", type=str, default='98')

    args = parser.parse_args()
    subloc = SublocProteinPipeline(args.fasta, args.denovo, args.species, args.go_list, args.out, args.go_identity)
    subloc.run()