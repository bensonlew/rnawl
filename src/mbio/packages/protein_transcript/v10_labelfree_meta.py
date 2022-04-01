#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/26 17:51
@file    : labelfree_v10_meta.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/labelfree_v10')
from v10_labelfree import V10Labelfree
sys.path.insert(0, '/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
from mako.template import Template
import os
import argparse


class V10LabelfreeMeta(V10Labelfree):
    def __init__(self, config, extract_nr, idmappingDB):
        super(V10LabelfreeMeta, self).__init__(config)
        self.idmappingDB = os.path.abspath(idmappingDB)
        self.extract_nr = True
        if extract_nr.lower() == 'no':
            self.extract_nr = False

    def add_commands(self):
        super(V10LabelfreeMeta, self).add_commands()
        params = dict(run_wd=self.output)
        params_ = dict(run_wd=self.output, qsub=0)
        self.taxo = self.add_command(name='taxo', params=params_)
        self.taxo_diff = self.add_command(name='taxo_diff', params=params_)
        self.pca = self.add_command(name='pca', params=params)

    def how_run(self):
        self.after_one(self.annot_diomand, self.run_go_result_nr)
        self.after_one(self.go_result_nr, self.run_go_result)
        self.after_one(self.annot_diomand, self.run_kegg_result)
        self.after_one(self.annot_diomand, self.run_cog_result)
        # self.after_one(self.annot_diomand, self.run_taxo)
        # 为了避免其跑错目录的问题
        self.after_one(self.go_result_nr, self.run_taxo)
        self.after_one(self.taxo, self.run_pca)
        if self.run_ppi_:
            self.after_some([self.annot_diomand, self.diff], self.run_ppi)
        self.after_one(self.diff, self.run_cluster)
        self.after_one(self.cluster, self.run_venn)
        self.after_some([self.diff, self.annot_enrich], self.run_ipath)
        self.after_some([self.diff, self.kegg_result], self.run_ipath_picture)
        self.after_some([self.diff, self.taxo], self.run_taxo_diff)
        self.after_some([self.go_result, self.kegg_result, self.cog_result, self.diff], self.run_annot_enrich)
        col_list = [self.annot_enrich, self.qc, self.venn, self.ipath, self.taxo_diff]
        if self.run_ppi_:
            col_list.append(self.ppi)
        self.after_some(col_list, self.run_results_collect)
        self.after_one(self.results_collect, self.end)

    def run_taxo(self):
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/taxo/taxonomy_protein.py '
        cmd += '-acc_list ' + os.path.join(self.raw_data, 'exp.list')
        cmd += ' -exp_txt ' + os.path.join(self.raw_data, 'exp.txt')
        cmd += ' -sample_config ' + os.path.join(self.raw_data, 'group.txt')
        cmd += ' -idmappingDB ' + self.idmappingDB
        cmd += ' -out ' + os.path.join(self.taxo.work_dir, 'taxonomy_result')
        if self.extract_nr:
            cmd += ' -extract_nr ' + 'yes'
            cmd += ' -fasta_vs_nr ' + os.path.join(self.go_result_nr.work_dir, 'exp.fasta_vs_nr.fast.xls')
        else:
            cmd += ' -extract_nr ' + 'no'
        params = dict(
            cmd=cmd,
            node=1,
            memory=3
        )
        self.taxo.set_params(params)
        self.taxo.run()

    def run_taxo_diff(self):
        cmd = 'python /mnt/ilustre/users/yitong.feng/scripts/taxo/taxonomy_protein_diff.py '
        cmd += ' -group_file ' + os.path.join(self.raw_data, 'group.txt')
        cmd += ' -tax_file ' + os.path.join(self.taxo.work_dir, 'get_profile', 'nr.tax.xls')
        cmd += ' -diff_path ' + os.path.join(self.diff.work_dir, 'diff_result_for_other_analyse')
        cmd += ' -out ' + os.path.join(self.taxo_diff.work_dir, 'taxo_diff_result')
        if self.extract_nr:
            cmd += ' -exp_txt ' + os.path.join(self.taxo.work_dir, 'get_profile', 'dealed_exp.txt')
        else:
            cmd += ' -exp_txt ' + os.path.join(self.raw_data, 'exp.txt')
        cmd += ' -acc2acc ' + os.path.join(self.taxo.work_dir, 'get_accs', 'accs.list_relation')
        params = dict(
            cmd=cmd,
            node=1,
            memory=3
        )
        self.taxo_diff.set_params(params)
        self.taxo_diff.run()

    def run_pca(self):
        cmd = '{bin}/plot-pca.pl -i {profile} -o {out} -w 8 -h 8 -m {group} -g group'.format(
            bin=self.bin,
            profile=os.path.join(self.taxo.work_dir, 'get_profile', 'nr.profile.xls'),
            out=self.pca.work_dir,
            group=os.path.join(self.raw_data, 'group.txt')
        )
        params = dict(
            cmd=cmd,
            node=2,
            memory=6
        )
        self.pca.set_params(params)
        self.pca.run()

    def run_results_collect(self):
        des_file = os.path.join(self.raw_data, 'description.txt')
        if not os.path.exists(des_file):
            acc_lists = list()
            with open(os.path.join(self.annot_diomand.work_dir, 'diomand_results', 'exp.fasta_vs_nr.fast.xls'), 'r') as \
                fr, open(des_file, 'w') as fw:
                fw.write('Accession\tDescription\n')
                for line in fr:
                    if not line.strip():
                        continue
                    line = line.strip().split('\t')
                    acc, des = line[0], line[-1].split('>')[0]
                    if acc not in acc_lists:
                        fw.write(acc+'\t'+des+'\n')
                        acc_lists.append(acc)

        if self.single != 'T':
            MJ = 'MJ_PM_LABELFREE_nosingle_autoreporter_meta.py'
        else:
            MJ = 'MJ_PM_LABELFREE_nosingle_autoreporter_meta.py'
        cmd = Template(filename='/mnt/ilustre/users/yitong.feng/scripts/labelfree_v10/result_collecter_labelfree.bash')
        cmd = cmd.render(outdir=self.output,
                         cog=self.cog,
                         rawdata=self.raw_data,
                         qc_dir=self.qc.work_dir,
                         go_result=self.go_result.work_dir,
                         kegg_result=self.kegg_result.work_dir,
                         cog_result=self.cog_result.work_dir,
                         diff_result=os.path.join(self.diff.work_dir, 'diff_result'),
                         diff_result_=os.path.join(self.diff.work_dir, 'diff_result_for_other_analyse'),
                         ipath_result=self.ipath.work_dir,
                         ipath_picture=os.path.join(self.ipath_picture.work_dir, 'ipath_pictures'),
                         venn_result=self.venn.work_dir,
                         cluster_result=os.path.join(self.cluster.work_dir, 'cluster_result'),
                         annot_enrich_result=os.path.join(self.annot_enrich.work_dir, 'diff_annot_enrich_result'),
                         bin=self.bin,
                         up=self.up,
                         dup=self.dup,
                         MJ=MJ,
                         )
        cmd = cmd.split('\n\n\n\n')
        if self.run_ppi_:
            ppi_cmd = "mkdir -p 3.DiffExpAnalysis/3.7.ppi\n cp %s/* 3.DiffExpAnalysis/3.7.ppi" % os.path.join(self.ppi.work_dir, 'ppi_results')
            cmd[1] += ppi_cmd
        taxo_cmd = "mkdir -p 2.Annotation/{2.4.TAXONOMY,2.5.PCA/}\n"
        taxo_cmd += "mkdir -p 3.DiffExpAnalysis/3.8.Taxonomy\n"
        taxo_cmd += "cp -r %s/* 2.Annotation/2.4.TAXONOMY/\n"%os.path.join(self.taxo.work_dir, 'taxonomy_result')
        taxo_cmd += "cp -r %s/* 3.DiffExpAnalysis/3.8.Taxonomy/\n"%os.path.join(self.taxo_diff.work_dir, 'taxo_diff_result')
        taxo_cmd += "cp %s/{*pdf,*xls} 2.Annotation/2.5.PCA/\n"%self.pca.work_dir
        cmd[1] += '\n' + taxo_cmd
        cmd = '\n'.join(cmd)
        params = dict(
            cmd=cmd,
            node=3,
            memory=12
        )

        self.results_collect.set_params(params)
        self.results_collect.run()

    def run(self):
        self.add_commands()
        self.how_run()
        self.run_annot_diomand()
        self.run_diff()
        self.run_qc()
        self.fire()


if __name__ == '__main__':
    # __init__(self, config, extract_nr, idmappingDB):
    parser = argparse.ArgumentParser(description="The workflow to run labelfree and can run taxonomy")
    parser.add_argument("-config", type=str, required=True, help="the same as labelfree config")
    parser.add_argument("-extract_nr", type=str, default='no',
                        help='if your accession is not normal, you should choose yes')
    parser.add_argument("-idmappingDB", type=str, default='/mnt/ilustre/users//bing.yang/DB/GO/idmapping.tb',
                        help='the idmapping tab')

    args = parser.parse_args()
    LABELFREE_taxo = V10LabelfreeMeta(args.config, args.extract_nr, args.idmappingDB)
    LABELFREE_taxo.run()