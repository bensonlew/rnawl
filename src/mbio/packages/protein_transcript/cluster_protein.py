# -*- coding: utf-8 -*-
# author fengyitong 2019-02

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
import os
import argparse
import glob


class ClusterProtein(Controller):
    def __init__(self, exp, group_file, diff_path, clt, mr, ml, out):
        super(ClusterProtein, self).__init__()
        self.exp = os.path.abspath(exp)
        if not os.path.exists(self.exp):
            super(ClusterProtein, self).end(normal=1, out='你传入的表达量文件不存在')
        self.group_file = os.path.abspath(group_file)
        if not os.path.exists(self.group_file):
            super(ClusterProtein, self).end(normal=1, out='你传入的分组文件不存在')
        self.diff_path = os.path.abspath(diff_path)
        if not os.path.exists(self.diff_path):
            super(ClusterProtein, self).end(normal=1, out='你传入的diff结果路径不存在')
        self.DElists = glob.glob(self.diff_path + '/*.DE.list')
        if not self.DElists:
            super(ClusterProtein, self).end(normal=1, out='你传入的diff结果路径下没有DE文件')
        self.clt = clt
        self.mr = mr
        self.ml = ml
        self.out = os.path.abspath(out)
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(ClusterProtein, self).end(normal=1, out='输出文件夹路径不对')
        tmp = os.path.join(self.out, 'group')
        with open(self.group_file, 'r') as gr:
            header = gr.readline()
            if u'#' in header:
                group_info = gr.read()
                if header.startswith('#'):
                    tmp_ = ''
                    for line in group_info.split('\n'):
                        if line:
                            line = line.strip().split('\t')
                            if '' not in line:
                                tmp_ += line[1] + '\t' + line[0] + '\n'
                    group_info = tmp_
                with open(tmp, 'w') as gw:
                    gw.write(group_info)
                self.group_file = tmp
        self.diff_list = self.cat_diff()
        self.group_cluster = list()
        self.sample_cluster = list()

    def cat_diff(self):
        de_file = os.path.join(self.out, 'All.DE.list')
        dp_list = list()
        for de in self.DElists:
            with open(de, 'r') as dr:
                dp_list += dr.read().strip().split('\n')
        while '' in dp_list:
            dp_list.remove('')
        with open(de_file, 'w') as dw:
            dw.write("Accession\n")
            dw.write('\n'.join(list(set(dp_list))))
        return de_file

    def add_commands(self):
        self.get_sample_exp = self.add_command(name='sample_exp')
        self.get_group_exp = self.add_command(name='group_exp')
        for n in [5, 10, 20, 40, 80]:
            cluster_sample = self.add_command(name='cluster_sample_' + str(n))
            cluster_group = self.add_command(name='cluster_group_' + str(n))
            self.sample_cluster.append(cluster_sample)
            self.group_cluster.append(cluster_group)

    def how_run(self):
        self.after_one(self.get_sample_exp, self.run_sample_cluster)
        self.after_one(self.get_group_exp, self.run_group_cluster)
        self.after_some(self.sample_cluster + self.group_cluster, self.set_output)

    def end(self):
        super(ClusterProtein, self).end(out='your program finished')

    def run_get_sample_exp(self):
        cmd = '/mnt/ilustre/users/ting.kuang/ITRAQ/bin/get_exp_from_list.pl %s %s %s/sample.diffexp.txt'%(self.diff_list, self.exp, self.out)
        params = dict(
            cmd=cmd,
            node=1,
            memory=5
        )
        self.get_sample_exp.set_params(params)
        self.get_sample_exp.run()

    def run_get_group_exp(self):
        cmd = '/mnt/ilustre/users/ting.kuang/ITRAQ/bin/get_groupmean.py -e %s -s %s -o %s/group.exp.txt'%(self.exp, self.group_file, self.out) + '\n'
        cmd += '/mnt/ilustre/users/ting.kuang/ITRAQ/bin/get_exp_from_list.pl %s %s/group.exp.txt %s/group.diffexp.txt'%(self.diff_list, self.out, self.out)
        params = dict(
            cmd=cmd,
            node=1,
            memory=5
        )
        self.get_group_exp.set_params(params)
        self.get_group_exp.run()

    def run_sample_cluster(self):
        for cluster in self.sample_cluster:
            n = cluster.name.split('_')[-1]
            cmd = '/mnt/ilustre/users/ting.kuang/ITRAQ/bin/plot_heatmap_trendline.pl -i %s/sample.diffexp.txt -clt %s -n %s -o %s.Heatmap -mr %s -ml %s'%(self.out, self.clt, n, n, self.mr, self.ml)
            params = dict(
                cmd=cmd,
                node=1,
                memory=2
            )
            cluster.set_params(params)
            cluster.run()

    def run_group_cluster(self):
        for cluster in self.group_cluster:
            n = cluster.name.split('_')[-1]
            cmd = '/mnt/ilustre/users/ting.kuang/ITRAQ/bin/plot_heatmap_trendline.pl -i %s/group.diffexp.txt -clt %s -n %s -o %s.Heatmap -mr %s -ml %s'%(self.out, self.clt, n, n, self.mr, self.ml)
            params = dict(
                cmd=cmd,
                node=1,
                memory=2
            )
            cluster.set_params(params)
            cluster.run()

    def set_output(self):
        try:
            os.mkdir(os.path.join(self.out, 'All_samples'))
            os.mkdir(os.path.join(self.out, 'Grouped'))
        except:
            pass
        for cluster in self.sample_cluster:
            for file in os.listdir(cluster.work_dir):
                if file.endswith('.pdf'):
                    try:
                        os.link(os.path.join(cluster.work_dir, file), os.path.join(self.out, 'All_samples', file))
                    except:
                        pass
                if os.path.isdir(os.path.join(cluster.work_dir, file)):
                    # print('cp -r %s %s' %(os.path.join(cluster.work_dir, file), os.path.join(self.out, 'All_samples')))
                    try:
                        os.system('cp -r %s %s' %(os.path.join(cluster.work_dir, file), os.path.join(self.out, 'All_samples')))
                    except:
                        pass
        for cluster in self.group_cluster:
            for file in os.listdir(cluster.work_dir):
                if file.endswith('.pdf'):
                    try:
                        os.link(os.path.join(cluster.work_dir, file), os.path.join(self.out, 'Grouped', file))
                    except:
                        pass
                if os.path.isdir(os.path.join(cluster.work_dir, file)):
                    os.system('cp -r %s %s' %(os.path.join(cluster.work_dir, file), os.path.join(self.out, 'Grouped')))
        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_get_group_exp()
        self.run_get_sample_exp()
        self.fire()

if __name__ == '__main__':
    # def __init__(self, exp, group_file, diff_path, clt, mr, ml, out):
    parser = argparse.ArgumentParser(description="The script to run cluster analyse for protein")
    parser.add_argument("-exp", type=str, required=True, help="the express tab")
    parser.add_argument("-group_file", type=str, required=True, help="group file")
    parser.add_argument("-diff_path", type=str, required=True, help="the diff result path")
    parser.add_argument("-clt", type=str, default='both',help='cluster type:both,row,colum or none')
    parser.add_argument("-mr", type=str, default='8',help='right margin,default 8')
    parser.add_argument("-ml", type=str, default='8',help='low margin,default 8')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'cluster_results'))

    args = parser.parse_args()
    cluster = ClusterProtein(args.exp, args.group_file, args.diff_path, args.clt, args.mr, args.ml, args.out)
    cluster.run()