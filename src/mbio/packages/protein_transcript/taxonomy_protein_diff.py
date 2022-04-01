# -*- coding: utf-8 -*-
# author fengyitong 2019-02

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
import os
import argparse
import glob
import pandas as pd
from collections import defaultdict
import shutil


class TaxonomyProteinDiff(Controller):
    def __init__(self, diff_path, exp_txt, group_file, tax_file, out, acc2acc):
        super(TaxonomyProteinDiff, self).__init__()
        self.exp_txt = os.path.abspath(exp_txt)
        if not os.path.exists(self.exp_txt):
            super(TaxonomyProteinDiff, self).end(normal=1, out='你传入的表达量文件不存在')
        self.tax_file = os.path.abspath(tax_file)
        if not os.path.exists(self.tax_file):
            super(TaxonomyProteinDiff, self).end(normal=1, out='你传入的tax分析结果文件不存在')
        self.diff_path = os.path.abspath(diff_path)
        if not os.path.exists(self.diff_path):
            super(TaxonomyProteinDiff, self).end(normal=1, out='你传入的diff结果路径不存在')
        self.DElists = glob.glob(self.diff_path + '/*.DE.list')
        if not self.DElists:
            super(TaxonomyProteinDiff, self).end(normal=1, out='你传入的diff结果路径下没有DE文件')
        self.group_file = os.path.abspath(group_file)
        if not os.path.exists(self.group_file):
            super(TaxonomyProteinDiff, self).end(normal=1, out='你传入的分组文件不存在')

        col1 = list()
        col2 = list()
        with open(self.group_file, 'r') as sr:
            for line in sr:
                if line.strip():
                    _1, _2 = line.strip().split('\t')
                    col1.append(_1)
                    col2.append(_2)
        sample_list = col1
        group_list = col2
        if len(set(col1)) < len(set(col2)):
            sample_list,group_list = group_list,sample_list
        self.group2sample = defaultdict(list)
        for i in zip(sample_list,group_list):
            self.group2sample[i[1]].append(i[0])

        self.out = os.path.abspath(out)
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(TaxonomyProteinDiff, self).end(normal=1, out='输出文件夹路径不对')
        self.de_taxos = list()
    #     如果是idmapping的，需要把tax_file给转换一下
        if acc2acc and os.path.exists(acc2acc):
            with open(acc2acc) as ar,open('nr.tax.xls', 'w') as tw, open(tax_file) as tr:
                acc_dict = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in ar if line.strip()}
                # ar.seek(0,0)
                # acc_dict_ = {line.strip().split('\t')[1]:line.strip().split('\t')[0] for line in ar if line.strip()}
                tw.write(tr.readline())
                for line in tr:
                    if line.strip():
                        tmp = line.strip().split('\t')[0]
                        if tmp in acc_dict:
                            line = line.replace(tmp, acc_dict[tmp])
                            tw.write(line)
            self.tax_file = os.path.abspath('nr.tax.xls')
            with open(self.exp_txt) as er, open('exp.txt', 'w') as ew:
                ew.write(er.readline())
                for line in er:
                    if line.strip():
                        tmp = line.strip().split('\t')[0]
                        if tmp in acc_dict:
                            line = line.replace(tmp, acc_dict[tmp])
                            ew.write(line)
            self.exp_txt = os.path.abspath('exp.txt')

            # self.acc_dict = acc_dict_
        else:
            # self.acc_dict = dict()
            pass


    def add_commands(self):
        for de in self.DElists:
            name = os.path.basename(de).split('.DE.list')[0]
            taxo = self.add_command(name=name)
            self.de_taxos.append(taxo)

    def how_run(self):
        self.after_some(self.de_taxos, self.set_output)

    def end(self):
        super(TaxonomyProteinDiff, self).end(out='your program finished')

    def run_de_taxos(self):
        exp = pd.read_table(self.exp_txt, sep = '\t', index_col= 0)
        # taxo = pd.read_table(self.tax_file, sep = '\t', index_col= 0)
        for de_taxo in self.de_taxos:
            cmp = de_taxo.name
            other, control = cmp.split('_vs_')
            samples = self.group2sample[other] + self.group2sample[control]
            de_file = os.path.join(self.diff_path, cmp + '.DE.list')
            with open(de_file, 'r') as der:
                de_list = [line.strip() for line in der if line.strip()]
                # if self.acc_dict:
                #     de_list = [self.acc_dict[x] for x in de_list if x in self.acc_dict]
                    # print(de_list)
            de_exp = exp[samples]
            de_list = list(set(de_exp.index.tolist()) & set(de_list))
            de_exp = de_exp.loc[de_list, :]

            de_exp_out = os.path.join(de_taxo.work_dir, cmp+'_exp.txt')
            de_exp.to_csv(de_exp_out, sep='\t', index=True)
            cmd = '/mnt/ilustre/users/bingxu.liu/workspace/tabletools_select.pl -i {delist} -t {taxo} -n 1 |/mnt/ilustre/users/bingxu.liu/workspace/tabletools_add.pl -i {exp} -t - -n 1 |cut -f3- > {cmp}.nr.profile.xls.raw\n'.format(delist=de_file,taxo=self.tax_file,exp=de_exp_out, cmp=cmp)
            cmd += 'python {bin}/get_taxon_exp.py -e {cmp}.nr.profile.xls.raw -o {cmp}.nr.profile.xls\n'.format(bin='/mnt/ilustre/users/ting.kuang/ALL-SCRIPT',cmp=cmp)
            cmd += 'perl {bin}/taxLevel_profile.2.pl -i {cmp}.nr.profile.xls -o {cmp}.summary'.format(bin='/mnt/ilustre/users/ting.kuang/ALL-SCRIPT',cmp=cmp)

            params = dict(
                cmd=cmd,
                node=3,
                memory=10
            )
            de_taxo.set_params(params)
            de_taxo.run()

    def set_output(self):
        for de_taxo in self.de_taxos:
            de_summary = os.path.join(self.out, de_taxo.name + '.summary')
            if os.path.exists(de_summary):
                shutil.rmtree(de_summary)
            cmd = 'cp -r ' + os.path.join(de_taxo.work_dir, de_taxo.name + '.summary') + ' ' + de_summary
            os.system(cmd)
            for file in os.listdir(de_summary):
                if not file.endswith(('.xls', '.pdf')) or file == 'Rplots.pdf':
                    os.remove(os.path.join(de_summary, file))

        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_de_taxos()
        self.fire()

if __name__ == '__main__':
    # __init__(self, diff_path, exp_txt, group_file, tax_file, out):
    parser = argparse.ArgumentParser(description="The script to run de_taxos for protein")
    parser.add_argument("-exp_txt", type=str, required=True, help="the expression tab")
    parser.add_argument("-group_file", type=str, required=True, help="the group table")
    parser.add_argument("-tax_file", type=str, required=True, help="the tax result file")
    parser.add_argument("-diff_path", type=str, required=True, help="the diff result path")
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'de_taxo_results'))
    parser.add_argument("-acc2acc", type=str, default='')

    args = parser.parse_args()
    de_taxo_ = TaxonomyProteinDiff(args.diff_path, args.exp_txt, args.group_file, args.tax_file, args.out, args.acc2acc)
    de_taxo_.run()