# -*- coding: utf-8 -*-
# author fengyitong 2019-02

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')
from controller import Controller
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import os
import argparse
import glob
import shutil


class IpathProteinPipeline(Controller):
    def __init__(self, pathway_table, diff_path, out):
        super(IpathProteinPipeline, self).__init__()
        self.pathway_table = os.path.abspath(pathway_table)
        if not os.path.exists(self.pathway_table):
            super(IpathProteinPipeline, self).end(normal=1, out='你传入的pathway_table文件不存在')
        self.diff_path = os.path.abspath(diff_path)
        if not os.path.exists(self.diff_path):
            super(IpathProteinPipeline, self).end(normal=1, out='你传入的diff结果路径不存在')
        self.DElists = glob.glob(self.diff_path + '/*.DE.list')
        if not self.DElists:
            super(IpathProteinPipeline, self).end(normal=1, out='你传入的diff结果路径下没有DE文件')
        self.out = os.path.abspath(out)
        if os.path.exists(self.out):
            shutil.rmtree(self.out)
        try:
            os.mkdir(self.out)
        except:
            super(IpathProteinPipeline, self).end(normal=1, out='输出文件夹路径不对')
        self.ipath_list = list()
        # self.script = '/mnt/ilustre/centos7users/yitong.feng/script/ipath/ipath_protein.py'
        self.script = '/mnt/ilustre/centos7users/yitong.feng/script/ipath/ipath3_protein.py'

    def add_commands(self):
        for de in self.DElists:
            de = os.path.basename(de)
            name = de.split('.DE.list')[0]
            ipath = self.add_command(name=name)
            self.ipath_list.append(ipath)

    def how_run(self):
        self.after_some(self.ipath_list, self.set_output)

    def end(self):
        super(IpathProteinPipeline, self).end(out='your program finished')

    def run_ipath(self):
        for ipath in self.ipath_list:
            name = ipath.name
            de = self.diff_path + '/' + name + '.DE.list'
            up = self.diff_path + '/' + name + '.up.list'
            down = self.diff_path + '/' + name + '.down.list'
            cmd = """python ${script} ${de} ${pathway_table} ${name}.DE
python ${script} ${up} ${pathway_table} ${name}.up
python ${script} ${down} ${pathway_table} ${name}.down
python ${script} ${up},${down} ${pathway_table} ${name}
"""
            cmd = Template(cmd)
            cmd = cmd.render(script=self.script,
                               up=up,
                               down=down,
                               de=de,
                               pathway_table=self.pathway_table,
                               name=ipath.name,
                               )
            params = dict(
                cmd=cmd,
                node=2,
                memory=6
            )
            ipath.set_params(params)
            ipath.run()

    def set_output(self):
        for ipath in self.ipath_list:
            for file in os.listdir(ipath.work_dir):
                if file.endswith(('.pdf','.svg','.png')):
                    try:
                        os.link(os.path.join(ipath.work_dir, file), os.path.join(self.out, file))
                    except:
                        pass
        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_ipath()
        self.fire()

if __name__ == '__main__':
    # def __init__(self, pathway_table, diff_path, out):
    parser = argparse.ArgumentParser(description="The script to create ipath pictures for protein")
    parser.add_argument("-pathway_table", type=str, required=True, help="the pathway_table from kegg annot result")
    parser.add_argument("-diff_path", type=str, required=True, help="the diff result path")
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'ipath_results'))

    args = parser.parse_args()
    ipath = IpathProteinPipeline(args.pathway_table, args.diff_path, args.out)
    ipath.run()