#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import shutil
import os
import pandas as pd
from collections import defaultdict
import argparse
import glob


class DiffLabelfree(Controller):
    def __init__(self, exp, group_file, control_file, method, alternative, p_correct, up, down, pvalue, out, cutoffs,average):
        super(DiffLabelfree, self).__init__()
        self.exp = os.path.abspath(exp)
        if not os.path.exists(self.exp):
            super(DiffLabelfree, self).end(normal=1, out='你传入的表达量文件不存在')
        # 处理exp文件，删除其所有列都为空的行
        exp_pd = pd.read_csv(self.exp,sep='\t',index_col=0)
        e = exp_pd.dropna(how='all', subset=exp_pd.columns.tolist())
        e.to_csv(self.exp, sep='\t', header=True, index=True)
        if not os.path.exists(group_file):
            super(DiffLabelfree, self).end(normal=1, out='你传入的分组文件不存在')
        if not os.path.exists(control_file):
            super(DiffLabelfree, self).end(normal=1, out='你传入的对比组文件不存在')
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
                    super(DiffLabelfree, self).end(normal=1, out='分组文件中存在表达量文件里没有的样本%s'%sample)
        self.groups = df_group['group'].tolist()
        for i in self.control2other:
            if i[0] not in self.groups or i[1] not in self.groups:
                super(DiffLabelfree, self).end(normal=1, out='分组文件与对照组文件不匹配')
        self.diff_commands = list()
        self.alternative = alternative
        self.method = method
        self.p_correct = p_correct
        self.up = up
        self.down = down
        self.pvalue = pvalue
        self.out = os.path.abspath(out)
        if os.path.exists(self.out):
            shutil.rmtree(self.out)
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(DiffLabelfree, self).end(normal=1, out='输出文件夹路径不对')
        self.g2cutoff = dict()
        if cutoffs != 'helpyourself':
            cutoffs = cutoffs.split(',')
            if len(self.groups) != len(cutoffs):
                super(DiffLabelfree, self).end(normal=1, out='指定cutoff的数量跟分组的数量不匹配')
            else:
                self.g2cutoff = {x[1]:x[2] for x in zip(self.groups, cutoffs)}
        else:
            self.g2cutoff = {g:float(len(self.group2sample[g]))/2 for g in self.groups}
        self.average = average.lower()


    def add_commands(self):
        for cmp in self.control2other:
            name = cmp[1] + '_vs_' + cmp[0]
            diff_c = self.add_command(name=name)
            self.diff_commands.append(diff_c)

    def how_run(self):
        self.after_some(self.diff_commands, self.set_output)

    def end(self):
        super(DiffLabelfree, self).end(out='your program finished')

    def run_diff(self):
        for diff in self.diff_commands:
            other_name, control_name = diff.name.split('_vs_')
            other = ','.join([str(self.samples.index(x)+1) for x in self.group2sample[other_name]])
            control = ','.join([str(self.samples.index(x)+1) for x in self.group2sample[control_name]])
            out = os.path.join(diff.work_dir, diff.name + '.diff.exp.xls')
            sympol = diff.name.replace('_vs_', '/')
            # r_cmd = open('/mnt/ilustre/users/yitong.feng/scripts/diff/diff_group.r', 'r').read()
            r_cmd = Template(filename='/mnt/ilustre/users/yitong.feng/scripts/diff/diff_labelfree.r',input_encoding='utf-8',output_encoding='utf-8')
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
                                 cutoff_a=self.g2cutoff[control_name],
                                 cutoff_b=self.g2cutoff[other_name],
                                 up=float(self.up),
                                 down=float(self.down),
                                 p_cutoff=float(self.pvalue),
                                 average=self.average
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
        dl = self.out + "_for_other_analyse"
        if os.path.exists(dl):
            shutil.rmtree(dl)
        os.mkdir(dl)
        for diff in self.diff_commands:
            for file in os.listdir(diff.work_dir):
                if file.endswith('.xls') or file.endswith('.list') or file.endswith('.pdf'):
                    target = os.path.join(self.out, file)
                    if os.path.exists(target):
                        os.remove(target)
                    # os.link(os.path.join(diff.work_dir, file), target)
                    cmd = 'cp ' + os.path.join(diff.work_dir, file) + ' ' + target
                    os.system(cmd)
            with open(os.path.join(diff.work_dir, diff.name+'.diff.exp.xls'), 'r') as diff_r,\
                open(os.path.join(diff.work_dir, diff.name+'.10.xls'), 'r') as _10_r,\
                open(os.path.join(diff.work_dir, diff.name+'.01.xls'), 'r') as _01_r, \
                open(os.path.join(diff.work_dir, diff.name+'.00.xls'), 'r') as _00_r, \
                open(os.path.join(dl, diff.name+'.diff.exp.xls'), 'w') as diff_w,\
                open(os.path.join(dl, diff.name+'.DE.list'), 'w') as de_w, \
                open(os.path.join(dl, diff.name+'.up.list'), 'w') as up_w, \
                open(os.path.join(dl, diff.name+'.down.list'), 'w') as down_w:
                for line in diff_r:
                    if line.strip():
                        diff_w.write(line.strip() + '\n')
                        tmp = line.strip().split('\t')
                        if 'up' in tmp and 'yes' in tmp:
                            de_w.write(tmp[0] + '\n')
                            up_w.write(tmp[0] + '\n')
                        if 'down' in tmp and 'yes' in tmp:
                            de_w.write(tmp[0] + '\n')
                            down_w.write(tmp[0] + '\n')
                _ = _01_r.readline()
                for line in _01_r:
                    if line.strip():
                        diff_w.write(line.strip() + '\n')
                        tmp = line.strip().split('\t')
                        de_w.write(tmp[0] + '\n')
                        up_w.write(tmp[0] + '\n')

                _ = _10_r.readline()
                for line in _10_r:
                    if line.strip():
                        diff_w.write(line.strip() + '\n')
                        tmp = line.strip().split('\t')
                        de_w.write(tmp[0] + '\n')
                        down_w.write(tmp[0] + '\n')

        def get_line_num(file):
            line_num = 0
            with open(file, 'r') as file_r:
                for line in file_r.readlines():
                    if len(line) and 'ccession' not in line:
                        line_num += 1
            return line_num

        file_list = glob.glob(self.out + '/*.DE.list')
        file_name = [file.strip().split('.DE.list')[0] for file in file_list]
        with open(os.path.join(self.out, 'all_diff_up_down.xls'), 'w')as file_each_a:
            file_each_a.write("name(A_vs_B)\tall\tboth\tdiff_num\tup_num\tdown_num\tonly A\tonly B\tneither\n")
            for each_name in file_name:
                both = get_line_num('%s.diff.exp.xls' % each_name)
                diff = get_line_num('%s.DE.list' % each_name)
                up = get_line_num('%s.up.list' % each_name)
                down = get_line_num('%s.down.list' % each_name)
                only_A = get_line_num('%s.10.xls' % each_name)
                only_B = get_line_num('%s.01.xls' % each_name)
                neither = get_line_num('%s.00.xls' % each_name)
                either = only_A + only_B + neither
                all = '%s=%s+%s' % (both + either, both, either)
                num_list = map(str, [all, both, diff, up, down, only_A, only_B, neither])
                file_each_a.write(os.path.basename(each_name) + '\t' + '\t'.join(num_list) + '\n')
        # 生成venn分析前体文件
        venn_file = os.path.join(dl, 'venn_pre')
        if os.path.exists(venn_file):
            shutil.rmtree(venn_file)
        os.mkdir(venn_file)
        with open(self.exp, 'r') as e_r:
            header = e_r.readline().strip().split('\t')
            for line in e_r:
                if line.strip():
                    line = line.strip().split('\t')
                    for s in header[1:]:
                        n = header.index(s)
                        try:
                            e = float(line[n])
                            if e > 0:
                                with open(os.path.join(venn_file, '%s.venn.list'%s), 'a') as fw:
                                    fw.write(line[0] + '\n')
                        except:
                            pass
        # 把结果调成跟线下的一致
        exp_pd = pd.read_csv(self.exp, index_col=0,sep='\t')
        for diff in self.diff_commands:
            other_name, control_name = diff.name.split('_vs_')
            other = self.group2sample[other_name]
            control = self.group2sample[control_name]
            with open(os.path.join(diff.work_dir, diff.name+'.10.xls'), 'r') as _10_r, \
                    open(os.path.join(diff.work_dir, diff.name + '.00.xls'), 'r') as _00_r,\
                    open(os.path.join(diff.work_dir, diff.name + '.01.xls'), 'r') as _01_r:
                _ = _10_r.readline()
                _ = _01_r.readline()
                _ = _00_r.readline()
                _10_list = [line.split('\t')[0] for line in _10_r if line.strip()]
                _01_list = [line.split('\t')[0] for line in _01_r if line.strip()]
                _00_list = [line.split('\t')[0] for line in _00_r if line.strip()]
                # print(_10_list)
                pd_10 = exp_pd.loc[_01_list, other+control]
                pd_10.to_csv(os.path.join(self.out, diff.name+'.10.xls'),sep='\t', header=True,
                               index=True)
                pd_01 = exp_pd.loc[_10_list, other + control]
                pd_01.to_csv(os.path.join(self.out, diff.name + '.01.xls'), sep='\t', header=True,
                             index=True)
                pd_00 = exp_pd.loc[_00_list, other + control]
                pd_00.to_csv(os.path.join(self.out, diff.name + '.00.xls'), sep='\t', header=True,
                             index=True)

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
    parser.add_argument("-cutoffs", type=str, default='helpyourself')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'diff_results'))
    parser.add_argument("-average", type=str, default='false', help='whether to average the zero')

    args = parser.parse_args()
    diff = DiffLabelfree(args.exp, args.group_file, args.control_file, args.method, args.alternative, args.p_correct, args.up,
                         args.down, args.pvalue, args.out, args.cutoffs, args.average)
    diff.run()