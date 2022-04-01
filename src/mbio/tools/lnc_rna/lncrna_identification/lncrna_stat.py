#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/5 16:18
@file    : lncrna_stat.py
"""
import csv
import json
import os
import unittest
from collections import defaultdict
from itertools import chain

import pandas as pd

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class LncrnaStatAgent(Agent):
    def __init__(self, parent):
        super(LncrnaStatAgent, self).__init__(parent)
        options = [
            {'name': 'novel_lncrna_detail', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'known_lncrna_detail', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'gene_type', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'new_gene_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'new_gene_list', 'type': 'infile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps("lncrna_stat")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.lncrna_stat.start()
        self.step.update()

    def step_end(self):
        self.step.lncrna_stat.finish()
        self.step.update()

    def check_options(self):
        for name in ('novel_lncrna_detail', 'known_lncrna_detail', 'exp_matrix', 'gene_type'):
            if not self.option(name).is_set:
                raise OptionError('缺少 %s 文件' % name)
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(LncrnaStatAgent, self).end()


class LncrnaStatTool(Tool):
    def new_genes(self):
        new_gene_set = set()
        with self.option('new_gene_gtf') as gtf_handler:
            for line_list, _ in gtf_handler:
                attr_dict = line_list[8]
                gene_id = attr_dict['gene_id']
                new_gene_set.add(gene_id)
        return new_gene_set

    def parser_detail_info(self, reader, file_type='novel'):
        # known lncrna_id       gene_id chromosome      ...       type    database
        # novel transcript_id   gene_id gene_name       ...      type
        trans_name = 'transcript_id' if file_type == 'novel' else 'lncrna_id'
        res_dict = {}
        for dic in reader:
            res_dict[dic[trans_name]] = [dic['gene_id'], dic['type']]
        return res_dict

    def gene_types(self):
        res_dict = defaultdict(set)
        for gene, _, rna_type, k_type in self.option('gene_type').csv_reader(remove_header=False):
            res_dict[(rna_type, k_type)].add(gene)
        return res_dict

    def statistic(self, known_dic, novel_dic, new_genes_set, gene_types):
        if not isinstance(new_genes_set, set):
            new_genes_set = set(new_genes_set)  # 废弃参数

        k_t_set = {i for i in known_dic}
        n_t_set = {i for i in novel_dic}

        exp_df = self.option('exp_matrix').dataframe(header=0, index_col=0)
        columns = exp_df.columns

        m_k_g_set = gene_types[('mRNA', 'known')]  # 已知mRNA基因
        m_n_g_set = gene_types[('mRNA', 'novel')]  # 新mRNA基因
        l_k_g_set = gene_types[('lncRNA', 'known')]  # 已知lncRNA基因
        l_n_g_set = gene_types[('lncRNA', 'novel')]  # 新lncRNA基因

        sample_stat_file = os.path.join(self.output_dir, 'lncrna_stat_in_sample.xls')
        category_stat_file = os.path.join(self.output_dir, 'lncrna_stat_in_category.xls')
        with open(sample_stat_file, 'w') as s_handler:
            fields = ['sample_name', 'known_num', 'new_num', 'total_num', 'type']
            line_demo = '\t'.join('{%s}' % i for i in fields) + '\n'
            s_handler.write('\t'.join(fields) + '\n')

            # s: sanple, k: known, t: transcript, g: gene
            for sample in columns:
                s_ser = exp_df[sample]
                trans_set = {i for i in s_ser[s_ser > 0].index}
                s_k_t_set = trans_set & k_t_set  # 已知lncRNA
                s_n_t_set = trans_set & n_t_set  # 新lncRNA
                s_handler.write(line_demo.format(sample_name=sample, known_num=len(s_k_t_set), new_num=len(s_n_t_set),
                                                 total_num=len(s_k_t_set) + len(s_n_t_set), type='LT'))

                all_genes = {known_dic[i][0] for i in s_k_t_set} | {novel_dic[i][0] for i in s_n_t_set}
                s_l_k_g_set = all_genes & l_k_g_set  # 已知lncRNA基因
                s_l_n_g_set = all_genes & l_n_g_set  # 新lncRNA基因
                s_handler.write(line_demo.format(sample_name=sample,
                                                 known_num=len(s_l_k_g_set), new_num=len(s_l_n_g_set),
                                                 total_num=len(s_l_k_g_set) + len(s_l_n_g_set), type='LG'))
                s_m_k_g_set = all_genes & m_k_g_set  # 已知mRNA基因
                s_m_n_g_set = all_genes & m_n_g_set  # 新mRNA基因
                s_handler.write(line_demo.format(sample_name=sample,
                                                 known_num=len(s_m_k_g_set), new_num=len(s_m_n_g_set),
                                                 total_num=len(s_m_k_g_set) + len(s_m_n_g_set), type='G'))

        s_stat = defaultdict(list)
        for tran_id, (gene, c_type) in chain(known_dic.items(), novel_dic.items()):
            s_stat[c_type].append((tran_id, gene))

        with  open(category_stat_file, 'w') as c_handler:
            fields = ['biotype', 'novel', 'known', 'total', 'type']
            line_demo = '\t'.join('{%s}' % i for i in fields) + '\n'
            c_handler.write('\t'.join(fields) + '\n')

            # c: category, m: mRNA, l: lncRNA, k: known, n: novel, t: transcript, g: gene
            for c_type, seq_list in s_stat.items():
                t_set, g_set = [set(i) for i in zip(*seq_list)]
                c_k_t_set = t_set & k_t_set  # 已知lncRNA
                c_n_t_set = t_set & n_t_set  # 新lncRNA
                c_handler.write(line_demo.format(biotype=c_type, novel=len(c_n_t_set), known=len(c_k_t_set),
                                                 total=len(c_k_t_set) + len(c_n_t_set), type='LT'))

                c_l_k_g_set = g_set & l_k_g_set  # 已知lncRNA基因
                c_l_n_g_set = g_set & l_n_g_set  # 新lncRNA基因
                c_handler.write(line_demo.format(biotype=c_type,
                                                 known=len(c_l_k_g_set), novel=len(c_l_n_g_set),
                                                 total=len(c_l_k_g_set) + len(c_l_n_g_set), type='LG'))

                c_m_k_g_set = g_set & m_k_g_set  # 已知mRNA基因
                c_m_n_g_set = g_set & m_n_g_set  # 新mRNA基因
                c_handler.write(line_demo.format(biotype=c_type,
                                                 known=len(c_m_k_g_set), novel=len(c_m_n_g_set),
                                                 total=len(c_m_k_g_set) + len(c_m_n_g_set), type='G'))

    def run(self):
        super(LncrnaStatTool, self).run()
        if self.option('new_gene_list').is_set:
            new_genes_set = set(map(str.strip, open(self.option('new_gene_list').path)))
        else:
            new_genes_set = self.new_genes()
        gene_types = self.gene_types()
        known_detail = self.parser_detail_info(self.option('known_lncrna_detail').csv_dict_reader(), file_type='known')
        novel_detail = self.parser_detail_info(self.option('novel_lncrna_detail').csv_dict_reader(), file_type='novel')
        self.statistic(known_detail, novel_detail, new_genes_set, gene_types)
        self.end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet

            data = {
                "id": "lncrna_stat_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_identification.lncrna_stat",
                "instant": False,
                "options": dict(
                    # {'name': 'novel_lncrna_detail', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
                    novel_lncrna_detail='/mnt/ilustre/users/sanger-dev/workspace/20190411/'
                                        'Single_new_lncrna_predict3512/NewLncrnaPredict/output/'
                                        'novel_lncrna_predict_detail.xls',
                    # {'name': 'known_lncrna_detail', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
                    known_lncrna_detail='/mnt/ilustre/users/sanger-dev/workspace/20190411/'
                                        'Single_new_lncrna_predict1150/KnownLncIdentify/output/known_lncrna_detail.xls',
                    # {'name': 'exp_matrix', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
                    exp_matrix='/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/transcript.tpm.matrix',
                    # {'name': 'new_gene_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
                    new_gene_gtf='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-cufflinks/'
                                 'output/NewTranscripts/new_genes.gtf',
                    gene_type=''
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
