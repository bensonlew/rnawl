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
import pandas as pd
from collections import defaultdict
import argparse


class DiffItraq(Controller):
    def __init__(self, exp, group_file, control_file, method, alternative, p_correct, up, down, pvalue, out):
        super(DiffItraq, self).__init__()
        self.exp = os.path.abspath(exp)
        if not os.path.exists(self.exp):
            super(DiffItraq, self).end(normal=1, out='你传入的表达量文件不存在')
        if not os.path.exists(group_file):
            super(DiffItraq, self).end(normal=1, out='你传入的分组文件不存在')
        if not os.path.exists(control_file):
            super(DiffItraq, self).end(normal=1, out='你传入的对比组文件不存在')
        with open(self.exp, 'r') as e_r:
            header = e_r.readline()
            self.samples = header.strip().split('\t')[1:]
        df_group = pd.read_table(group_file, dtype={0:str,1:str})
        df_control = pd.read_table(control_file, dtype={0:str,1:str})
        self.control2other = zip(df_control['#control'],df_control['other'])
        sample2group = zip(df_group['#sample'], df_group['group'])
        self.group2sample = defaultdict(list)
        for i in sample2group:
            self.group2sample[i[1]].append(i[0])
        for group in self.group2sample:
            for sample in self.group2sample[group]:
                if sample not in self.samples:
                    super(DiffItraq, self).end(normal=1, out='分组文件中存在表达量文件里没有的样本%s'%sample)
        self.groups = df_group['group'].tolist()
        for i in self.control2other:
            if i[0] not in self.groups or i[1] not in self.groups:
                super(DiffItraq, self).end(normal=1, out='分组文件与对照组文件不匹配')
        self.diff_commands = list()
        self.alternative = alternative
        self.method = method
        self.p_correct = p_correct
        self.up = up
        self.down = down
        self.pvalue = pvalue
        self.out = out
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(DiffItraq, self).end(normal=1, out='输出文件夹路径不对')

    def add_commands(self):
        for cmp in self.control2other:
            name = cmp[1] + '_vs_' + cmp[0]
            diff_c = self.add_command(name=name)
            self.diff_commands.append(diff_c)

    def how_run(self):
        self.after_some(self.diff_commands, self.set_output)

    def end(self):
        super(DiffItraq, self).end(out='your program finished')

    def run_diff(self):
        for diff in self.diff_commands:
            other_name, control_name = diff.name.split('_vs_')
            other = ','.join([str(self.samples.index(x)+1) for x in self.group2sample[other_name]])
            control = ','.join([str(self.samples.index(x)+1) for x in self.group2sample[control_name]])
            out = os.path.join(diff.work_dir, diff.name + '.diff.exp.xls')
            sympol = diff.name.replace('_vs_', '/')
            # r_cmd = open('/mnt/ilustre/users/yitong.feng/scripts/diff/diff_group.r', 'r').read()
            r_cmd = Template(filename='/mnt/ilustre/users/yitong.feng/scripts/diff/diff_group.r',input_encoding='utf-8',output_encoding='utf-8')
            # r_cmd = Template(r_cmd)
            r_cmd = r_cmd.render(exp=self.exp,
                                 control=control,
                                 other=other,
                                 method=self.method,
                                 alternative=self.alternative,
                                 p_correct=self.p_correct,
                                 other_name=other_name,
                                 control_name=control_name,
                                 sympol=sympol,
                                 outxls=out,
                                 up=float(self.up),
                                 down=float(self.down),
                                 p_cutoff=float(self.pvalue)
                                 )
            r_path = os.path.join(diff.work_dir, diff.name + 'diff.r')
            with open(r_path, 'w') as rw:
                rw.write(r_cmd)
            cmd = 'Rscript ' + r_path
            cmd_ = r"""less ${name}.diff.exp.xls | awk '{if($6<${pvalue} && ($4>${up} || $4<${down})){print $1}}' > ${name}.DE.list
less ${name}.diff.exp.xls | awk '{if($6<${pvalue} && $4>${up}){print $1}}' > ${name}.up.list
less ${name}.diff.exp.xls | awk '{if($6<${pvalue} && $4<${down}){print $1}}' > ${name}.down.list
less ${name}.diff.exp.xls | awk -F "\t" '{printf $1"\t"$5"\t"$6"\t"; if(NR==1){print "sig"}else if($6<0.05){if($4>${up}){printf "up"; if($6<0.01){print "-p-0.01"}else{print "-p-0.05"}}else if($4<${down}){printf "down"; if($6<0.01){print "-p-0.01"}else{print "-p-0.05"}}else{print "nosig"}}else{print "nosig"}}' > ${name}.volcano
less ${name}.diff.exp.xls | awk -F "\t" '{printf $1"\t"$2"\t"$3"\t"; if(NR==1){print "sig"}else if($6<0.05){if($4>${up}){printf "up"; if($6<0.01){print "-p-0.01"}else{print "-p-0.05"}}else if($4<${down}){printf "down"; if($6<0.01){print "-p-0.01"}else{print "-p-0.05"}}else{print "nosig"}}else{print "nosig"}}' > ${name}.scatter
/mnt/ilustre/users/ting.kuang/ITRAQ/bin/Highchart_for_ITRAQ.pl -yAxis_log -yAxis_min 0.0001 -yAxis_max 1 -type scatter -t ${name}.volcano -yAxis_reversed -scatter_series down-p-0.01,down-p-0.05,nosig,up-p-0.05,up-p-0.01 -width 700 -height 500 -color_the "'#2222FF','#22CCFF','#222222','#FFCC22','#FF2222'" -scatter_size 2,2,1,2,2 -scatter_symbol "'triangle-down','triangle-down','diamond','triangle','triangle'"
/mnt/ilustre/users/ting.kuang/ITRAQ/bin/Highchart_for_ITRAQ.pl -type scatter -t ${name}.scatter -xAxis_log -yAxis_log -scatter_series down-p-0.01,down-p-0.05,nosig,up-p-0.05,up-p-0.01 -width 700 -height 500 -color_the "'#2222FF','#22CCFF','#222222','#FFCC22','#FF2222'" -scatter_size 2,2,2,2,2 -scatter_symbol "'triangle-down','triangle-down','diamond','triangle','triangle'"
"""
            cmd_ = Template(cmd_)
            cmd_ = cmd_.render(name=diff.name,
                                 up=self.up,
                                 down=self.down,
                                 pvalue=self.pvalue,
                                 p_correct=self.p_correct,
                                 )
            cmd += '\n' + cmd_
            params = dict(
                cmd=cmd,
                node=2,
                memory=6
            )
            diff.set_params(params)
            diff.run()

    def set_output(self):
        for diff in self.diff_commands:
            for file in os.listdir(diff.work_dir):
                if file.endswith('.xls') or file.endswith('.list') or file.endswith('.pdf'):
                    target = os.path.join(self.out, file)
                    if os.path.exists(target):
                        os.remove(target)
                    os.link(os.path.join(diff.work_dir, file), target)
        os.chdir(self.out)
        os.system('/mnt/ilustre/users/ting.kuang/ITRAQ/bin/get_diff_up_down.py')
        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_diff()
        self.fire()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="The script to run diff analyse for itraq")
    parser.add_argument("-exp", type=str, required=True, help="the express tab")
    parser.add_argument("-group_file", type=str, required=True, help="group file")
    parser.add_argument("-control_file", type=str, required=True, help="control file")
    parser.add_argument("-method", type=str, default='student',help='"student", "welch", "wilcox", "chisq", "fisher"')
    parser.add_argument("-alternative", type=str, default='two.sided',help='"two.sided", "less", "greater"')
    parser.add_argument("-p_correct", type=str, default='bonferroni',help='"holm", "bonferroni", "BH", "BY"')
    parser.add_argument("-up", type=str, default='1.2')
    parser.add_argument("-down", type=str, default='0.83')
    parser.add_argument("-pvalue", type=str, default='0.05')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'diff_results'))

    args = parser.parse_args()
    diff = DiffItraq(args.exp, args.group_file, args.control_file, args.method, args.alternative, args.p_correct, args.up,
                         args.down, args.pvalue, args.out)
    diff.run()