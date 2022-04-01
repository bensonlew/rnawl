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


class PpiByBlast(Controller):
    def __init__(self, diff_path, db_file, blast_file, identity, out):
        super(PpiByBlast, self).__init__()
        self.db_file = os.path.abspath(db_file)
        if not os.path.exists(self.db_file):
            super(PpiByBlast, self).end(normal=1, out='你传入的db文件不存在')
        self.blast_file = os.path.abspath(blast_file)
        if not os.path.exists(self.blast_file):
            super(PpiByBlast, self).end(normal=1, out='你传入的blast结果文件不存在')
        self.diff_path = os.path.abspath(diff_path)
        if not os.path.exists(self.diff_path):
            super(PpiByBlast, self).end(normal=1, out='你传入的diff结果路径不存在')
        self.DElists = glob.glob(self.diff_path + '/*.DE.list')
        if not self.DElists:
            super(PpiByBlast, self).end(normal=1, out='你传入的diff结果路径下没有DE文件')
        self.identity = float(identity)
        self.out = os.path.abspath(out)
        if not os.path.exists(self.out):
            try:
                os.mkdir(self.out)
            except:
                super(PpiByBlast, self).end(normal=1, out='输出文件夹路径不对')
        self.ppis = list()

    def add_commands(self):
        for de in self.DElists:
            name = os.path.basename(de).split('.DE.list')[0]
            ppi = self.add_command(name=name+'_ppi')
            self.ppis.append(ppi)

    def how_run(self):
        self.after_some(self.ppis, self.set_output)

    def end(self):
        super(PpiByBlast, self).end(out='your program finished')

    def run_ppis(self):
        for ppi in self.ppis:
            cmp = ppi.name.split('_ppi')[0]
            de_list = os.path.join(self.diff_path, cmp + '.DE.list')
            dexls = os.path.join(self.diff_path, cmp + '.diff.exp.xls')
            info = os.path.join(ppi.work_dir, cmp+'.DE.info')
            with open(dexls, 'r') as de_r, open(info, 'w') as iw:
                for n, line in enumerate(de_r):
                    if line.strip():
                        tmp = line
                        line = line.strip().split('\t')
                        if n == 0:
                            iw.write(line[0] + '\t' + line[3] + '\n')
                        else:
                            if not 'no\t' in tmp and not 'no change\t' in tmp:
                                iw.write(line[0] + '\t' + line[3] + '\n')

            # cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/ppi/get_ppi_from_blast.py {} {} {} {} {}'.format(
            cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/ppi/get_ppi_from_blast_muilt.py {} {} {} {} {}'.format(
                self.blast_file, self.db_file, de_list, self.identity, cmp+'_ppi_result.xls'
            )
            # 加上直接出图的内容
            cmd += '\n' + 'python /mnt/ilustre/users/ting.kuang/ALL-SCRIPT/network_for_ppi.py -ppi_info {} -de_info {} -out {}'.format(cmp+'_ppi_result.xls', info, cmp)
            params = dict(
                cmd=cmd,
                node=9,
                memory=30
            )
            ppi.set_params(params)
            ppi.run()

    def set_output(self):
        for ppi in self.ppis:
            for file in os.listdir(ppi.work_dir):
                if file.endswith('.detail.xls'):
                    detail_xls = pd.read_table(os.path.join(ppi.work_dir, file), sep='\t')
                    del_list = [x for x in detail_xls.columns.tolist() if x.endswith('_transferred')]
                    detail_xls.drop(del_list, axis=1, inplace=True)
                    detail_xls = detail_xls.drop_duplicates()
                    detail_xls.to_csv(os.path.join(self.out, file), sep='\t', index=False)
                if u'.detail.' not in file and file.endswith('.xls'):
                    xls = pd.read_table(os.path.join(ppi.work_dir, file), sep='\t')
                    xls = xls.drop_duplicates()
                    xls.to_csv(os.path.join(self.out, file), sep='\t', index=False)
                if file.endswith(('.info', '.svg', '.pdf', '.png')):
                    try:
                        os.link(os.path.join(ppi.work_dir, file), os.path.join(self.out, file))
                    except:
                        pass
        self.end()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_ppis()
        self.fire()

if __name__ == '__main__':
    # def __init__(self, diff_path, db_file, blast_file, identity, out):
    parser = argparse.ArgumentParser(description="The script to run ppis for protein")
    parser.add_argument("-db_file", type=str, required=True, help="the ppi relation db")
    parser.add_argument("-blast_file", type=str, required=True, help="the blast result file, can be xls or xml")
    parser.add_argument("-diff_path", type=str, required=True, help="the diff result path")
    parser.add_argument("-identity", type=str, default='98',help='the identity cutoff')
    parser.add_argument("-out", type=str, default=os.path.join(os.getcwd(), 'ppi_results'))

    args = parser.parse_args()
    ppir = PpiByBlast(args.diff_path, args.db_file, args.blast_file, args.identity, args.out)
    ppir.run()