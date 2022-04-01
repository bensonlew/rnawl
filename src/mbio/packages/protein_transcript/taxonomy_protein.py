#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author fengyitong 2019-02

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import time
import os
from collections import defaultdict
import argparse
import pandas as pd


class TaxonomyProtein(Controller):
    def __init__(self, acc_list, exp_txt, fasta_vs_nr, sample_config, extract_nr, idmappingDB, out):
        super(TaxonomyProtein, self).__init__()
        self.acc_list = os.path.abspath(acc_list)
        if not os.path.exists(self.acc_list):
            super(TaxonomyProtein, self).end(normal=1, out='你传入的exp.list文件路径不正确')
        self.idmappingDB = os.path.abspath(idmappingDB)
        if not os.path.exists(self.idmappingDB):
            super(TaxonomyProtein, self).end(normal=1, out='你传入的idmappingDB文件路径不正确')
        self.extract_nr = True
        if extract_nr.lower() == 'no':
            self.extract_nr = False
        # with open(acc_list) as ar:
        #     for line in ar:
        #         if line.strip():
        #             if not line.strip().isalpha():
        #                 break
        #     else:
        #         self.extract_nr = False
        if fasta_vs_nr:
            self.fasta_vs_nr = os.path.abspath(fasta_vs_nr)
            if not os.path.exists(self.fasta_vs_nr):
                super(TaxonomyProtein, self).end(normal=1, out='你传入的nr比对结果文件路径不正确')
        if self.extract_nr and not fasta_vs_nr:
            super(TaxonomyProtein, self).end(normal=1, out='你的exp.list里的accession号不正规，需要跟nr的blast结果才能跑')
        self.exp_txt = os.path.abspath(exp_txt)
        if not os.path.exists(self.exp_txt):
            super(TaxonomyProtein, self).end(normal=1, out='你传入的exp_txt文件路径不正确')
        self.out = os.path.abspath(out)
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(TaxonomyProtein, self).end(normal=1, out='输出文件夹路径不对')
        if self.extract_nr:
            self.extract_nr_list()
        sample_config = os.path.abspath(sample_config)
        if not os.path.exists(sample_config):
            super(TaxonomyProtein, self).end(normal=1, out='你传入的分组文件路径不正确')
        self.sample_config = os.path.join(self.out, 'sample_config')
        col1 = list()
        col2 = list()
        with open(sample_config, 'r') as sr:
            for line in sr:
                if line.strip():
                    _1, _2 = line.strip().split('\t')
                    col1.append(_1)
                    col2.append(_2)
        sample_list = col1
        if len(set(col1)) < len(set(col2)):
            sample_list = col2
        with open(self.sample_config, 'w') as sw:
            for sample in sample_list:
                if not sample.startswith('#'):
                    sw.write(sample + '\t' + sample + '\n')

    def extract_nr_list(self):
        self.acc2nr = os.path.join(self.out, 'acc2nr.list')
        accs = list()
        with open(self.fasta_vs_nr, 'r') as nr, open(self.acc2nr, 'w') as nw:
            for line in nr:
                if line.strip():
                    line = line.strip().split('\t')
                    if line[0] not in accs:
                        nw.write(line[0] + '\t' + line[1] + '\n')
                        accs.append(line[0])

    def add_commands(self):
        if self.extract_nr:
            self.get_accs = self.add_command(name = 'get_accs')
        self.get_profile = self.add_command(name = 'get_profile')
        self.summary = self.add_command(name = 'summary')
        self.leveltree = self.add_command(name = 'leveltree')
        self.hierarchical = self.add_command(name = 'hierarchical')

    def how_run(self):
        if self.extract_nr:
            self.after_one(self.get_accs, self.run_get_profile)
        self.after_one(self.get_profile, self.run_summary)
        self.after_one(self.get_profile, self.run_leveltree)
        self.after_one(self.get_profile, self.run_hierarchical)
        self.after_some([self.summary, self.leveltree, self.hierarchical], self.set_output)

    def end(self):
        super(TaxonomyProtein, self).end(out='your program finished')

    def run_get_accs(self):
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/taxo/extract_accession_from_idmapping.py %s %s %s' %(
            self.acc2nr, self.idmappingDB, os.path.join(self.get_accs.work_dir, 'accs.list')
        )
        params = dict(
            cmd=cmd,
            node=9,
            memory=30
        )
        self.get_accs.set_params(params)
        self.get_accs.run()

    def deal_exp(self):
        with open(os.path.join(self.get_accs.work_dir, 'accs.list_relation'), 'r') as accr,\
            open(self.exp_txt, 'r') as expr,\
            open(os.path.join(self.get_profile.work_dir,'dealed_exp.txt'), 'w') as derp:
            acc2acc = {line.strip().split('\t')[1]:line.strip().split('\t')[0] for line in accr if line.strip()}
            derp.write(expr.readline())
            for line in expr:
                if line.strip():
                    tmp = line.strip().split('\t')
                    if tmp[0] in acc2acc:
                        derp.write(acc2acc[tmp[0]] + '\t' + '\t'.join(tmp[1:]) + '\n')
        self.exp_txt = os.path.join(self.get_profile.work_dir,'dealed_exp.txt')

    def run_get_profile(self):
        explist = self.acc_list
        if self.extract_nr:
            explist = os.path.join(self.get_accs.work_dir, 'accs.list')
            self.deal_exp()
        cmd = """python ${bin}/get_taxonomy_from_acclist.py -acc_file ${explist} -o nr.tax.xls
/mnt/ilustre/users/bingxu.liu/workspace/tabletools_add.pl -i ${exptxt} -t nr.tax.xls -n 1 |cut -f3- > nr.profile.xls.raw 
python ${bin}/get_taxon_exp.py -e nr.profile.xls.raw -o nr.profile.xls
"""
        cmd = Template(cmd)
        cmd = cmd.render(bin='/mnt/ilustre/users/ting.kuang/ALL-SCRIPT/',
                         explist=explist,
                         exptxt=self.exp_txt,
                         )
        params = dict(
            cmd=cmd,
            node=9,
            memory=30
        )
        self.get_profile.set_params(params)
        self.get_profile.run()

    def run_summary(self):
        profile = os.path.join(self.summary.work_dir, 'nr.profile.xls')
        try:
            os.link(os.path.join(self.get_profile.work_dir, 'nr.profile.xls'), profile)
        except:
            pass
        if not os.path.exists(profile):
            super(TaxonomyProtein, self).end(normal=1, out='get_profile那一步运行出错')
        cmd = 'perl /mnt/ilustre/users/ting.kuang/ALL-SCRIPT/taxLevel_profile.2.pl -i nr.profile.xls -o summary'
        params = dict(
            cmd=cmd,
            node=3,
            memory=10
        )
        self.summary.set_params(params)
        self.summary.run()

    def run_leveltree(self):
        profile = pd.read_table(os.path.join(self.get_profile.work_dir, 'nr.profile.xls'), sep='\t', index_col=0)
        profile = profile.fillna(0)
        samples = profile.columns.tolist()
        domain2tax = defaultdict(list)
        for i in profile.index:
            domain = i.split(';')[0].split('d__')[1]
            if domain != 'norank':
                domain2tax[domain].append(i)
        for domain in domain2tax:
            d_profile = profile.loc[domain2tax[domain], :]
            d_index = d_profile.index.tolist()
            d_profile['TaxID'] = ['Tax' + str(d_index.index(i)+1) for i in d_index]
            d_profile['Taxonomy'] = ['; '.join(ind.split(';')) for ind in d_profile.index]
            reorder = ['TaxID'] + samples + ['Taxonomy']
            d_profile = d_profile[reorder]
            d_profile.to_csv(os.path.join(self.leveltree.work_dir, domain+'_nr.profile.xls'), sep='\t', index=False)
        cmd = ''
        for domain in domain2tax:
            cmd += '/mnt/ilustre/users/ting.kuang/ALL-SCRIPT/level_tree.float.py -f %s -t 10 -g %s --hi 2000\n'%(os.path.join(self.leveltree.work_dir, domain+'_nr.profile.xls'), self.sample_config)
        params = dict(
            cmd=cmd,
            node=3,
            memory=10
        )
        self.leveltree.set_params(params)
        self.leveltree.run()

    def run_hierarchical(self):
        cmd = '/mnt/ilustre/users/ting.kuang/ALL-SCRIPT/krona.pl %s %s'%(os.path.join(self.get_profile.work_dir, 'nr.profile.xls'), self.hierarchical.work_dir)
        params = dict(
            cmd=cmd,
            node=3,
            memory=10
        )
        self.hierarchical.set_params(params)
        self.hierarchical.run()

    def set_output(self):
        summary = os.path.join(self.out, 'summary')
        if not os.path.exists(summary):
            os.mkdir(summary)
        for file in os.listdir(os.path.join(self.summary.work_dir, 'summary')):
            if file.endswith(('.xls', '.pdf')) and file != 'Rplots.pdf':
                source = os.path.join(self.summary.work_dir, 'summary', file)
                target = os.path.join(summary, file)
                if os.path.exists(target):
                    os.remove(target)
                os.link(source, target)

        level_tree = os.path.join(self.out, 'levelTree_multiSamples')
        if not os.path.exists(level_tree):
            os.mkdir(level_tree)
        for file in os.listdir(self.leveltree.work_dir):
            if file.endswith(('.xls', '.pdf', '.svg')):
                source = os.path.join(self.leveltree.work_dir, file)
                target = os.path.join(level_tree, file)
                if os.path.exists(target):
                    os.remove(target)
                os.link(source, target)

        hierarchical_pie = os.path.join(self.out, 'hierarchical_pie')
        if not os.path.exists(hierarchical_pie):
            os.mkdir(hierarchical_pie)
        for file in os.listdir(self.hierarchical.work_dir):
            if file.endswith(('.html', '.tax')):
                source = os.path.join(self.hierarchical.work_dir, file)
                target = os.path.join(hierarchical_pie, file)
                if os.path.exists(target):
                    os.remove(target)
                os.link(source, target)

        nr = os.path.join(self.out, 'nr.profile.xls')
        if os.path.exists(nr):
            os.remove(nr)
        os.link(os.path.join(self.get_profile.work_dir, 'nr.profile.xls'), nr)

        nr = os.path.join(self.out, 'nr.tax.xls')
        if os.path.exists(nr):
            os.remove(nr)
        os.link(os.path.join(self.get_profile.work_dir, 'nr.tax.xls'), nr)

        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        if self.extract_nr:
            self.run_get_accs()
        else:
            self.run_get_profile()
        self.fire()

if __name__ == '__main__':
    # __init__(self, acc_list, exp_txt, fasta_vs_nr, sample_config, extract_nr, idmappingDB, out):
    parser = argparse.ArgumentParser(description="The script to run taxonomy for protein")
    parser.add_argument("-acc_list", type=str, required=True, help="the accession list")
    parser.add_argument("-exp_txt", type=str, required=True, help="the expression tab")
    parser.add_argument("-sample_config", type=str, required=True, help="the group table")
    parser.add_argument("-fasta_vs_nr", type=str, default='',help='nr blast result')
    parser.add_argument("-extract_nr", type=str, default='no', help='if your accession is not normal, you should choose yes')
    parser.add_argument("-idmappingDB", type=str, default='/mnt/ilustre/users//bing.yang/DB/GO/idmapping.tb',help='the idmapping tab')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'taxonomy_result'))

    args = parser.parse_args()
    taxo = TaxonomyProtein(args.acc_list, args.exp_txt, args.fasta_vs_nr, args.sample_config, args.extract_nr, args.idmappingDB, args.out)
    taxo.run()