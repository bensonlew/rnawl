# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2019.09.15

from biocluster.module import Module
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
from biocluster.wsheet import Sheet
import os
import re
import shutil
import json
import glob
import pandas as pd
from biocluster.config import Config
import unittest

class TargetCistransModule(Module):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(TargetCistransModule, self).__init__(wsheet_object)
        options = [
            {"name": "novol", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {"name": "known", "type": "infile", "format": "lnc_rna.gtf"},  # 输入文件
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'annotation', 'format': 'lnc_rna.common'},
            {'name': 'exp_matrix_lnc', 'type': 'infile', 'format': 'lnc_rna.express_matrix'},
            {'name': 'exp_matrix_target', 'type': 'infile', 'format': 'lnc_rna.express_matrix'},
            {'name': 'pvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'qvalue_cutoff', 'type': 'float', 'default': 0.05},
            {'name': 'cor_cutoff', 'type': 'float', 'default': 0.8},
            {'name': 'corr_way', 'type': 'string', 'default': "spearmanr"},
            {'name': 'padjust_way', 'type': 'string', 'default': "fdr_bh"},
            {'type': 'string', 'name': 'geneset_id', 'default': 'All'},
            {'type': 'int', 'name': 'up_dis', 'default': 10},
            {'type': 'int', 'name': 'down_dis', 'default': 10},
            {'name': 'filter_lnc_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'filter_target_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {"name": "stat_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "string"},
            {"name": "target_cis_id", "type": "string"},
            {"name": "last_id_target", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "int"},
            {"name": "update_info", "type": "string"},
            {"name": "group_dict", "type": "string"},
            {"name": "group_id", "type": "string"},

        ]
        self.add_option(options)
        self.target_predict_known = self.add_tool("lnc_rna.lnc_target_cis")
        self.target_predict_novol = self.add_tool("lnc_rna.lnc_target_cis")
        self.merge_target = self.add_tool("lnc_rna.lnc_target_merge")
        self.expcorr_tool = self.add_tool("lnc_rna.expcorr_lnctarget")
        # self.target_predict.on('end', self.set_db)


    def get_seq_list(self, summary_file):
        a = pd.read_table(summary_file, header=0)
        return set(a.seq_id)

    def filter_exp(self):
        lnc_k_t2g = self.option("known").get_txpt_gene_dic()
        lnc_n_t2g = self.option("novol").get_txpt_gene_dic()
        m_t2g = self.option("mrna_gtf").get_txpt_gene_dic()


        lnc_diff_set = set()
        target_diff_set = set()
        if self.option("filter_lnc_file").is_set:
            lnc_diff_set = self.get_seq_list(self.option("filter_lnc_file").prop['path'])
        if self.option("filter_target_file").is_set:
            target_diff_set = self.get_seq_list(self.option("filter_target_file").prop['path'])

        if self.option("exp_matrix_lnc").is_set:
            # samples = self.option("use_samples").split("|")
            exp = pd.read_table(self.option("exp_matrix_lnc").prop['path'], header=0, index_col=0)
            if "rna_type" in exp.columns:
                exp = exp.drop(columns=["rna_type", "is_new"])
            # exp_max = exp.max(axis=1)
            # lnc_set = set(exp_max[exp_max>10].index)
            if len(lnc_diff_set) == 0:
                lnc_set = set(exp.index)
            else:
                lnc_set = set(exp.index).intersection(lnc_diff_set)

            exp_list = lnc_set
            gene_filter = list(set(lnc_k_t2g.keys() + lnc_n_t2g.keys()).intersection(set(exp_list)))
            exp_filter = exp.loc[gene_filter, :]
            exp_filter.to_csv(self.work_dir + "/lnc_matrix_filter.xls", sep="\t")
        else:
            pass

        if self.option("exp_matrix_target").is_set:
            # samples = self.option("use_samples").split("|")
            exp = pd.read_table(self.option("exp_matrix_target").prop['path'], header=0, index_col=0)
            if len(target_diff_set) == 0:
                exp_set = set(exp.index)
            else:
                exp_set = set(exp.index).intersection(target_diff_set)
            gene_filter = list(set(m_t2g.values()).intersection(exp_set))
            exp_filter = exp.loc[gene_filter, :]
            exp_filter.to_csv(self.work_dir + "/m_matrix_filter.xls", sep="\t")
        else:
            pass


    def run_exp_corr(self):
        options = dict(
            exp=self.work_dir + "/lnc_matrix_filter.xls",
            exp_t=self.work_dir + "/m_matrix_filter.xls",
            gt="G",
            anno=self.option('annotation').prop['path'],
            pvalue_cutoff=self.option('pvalue_cutoff'),
            qvalue_cutoff=self.option('qvalue_cutoff'),
            cor_cutoff=self.option('cor_cutoff'),
            corr_way=self.option('corr_way'),
            padjust_way=self.option('padjust_way'),
        )
        self.expcorr_tool.set_options(options)
        self.expcorr_tool.run()


    def run_known(self):
        options = {
            "lncrna_gtf" : self.option("known"),
            "mrna_gtf" : self.option("mrna_gtf"),
            "annotation" : self.option("annotation"),
            "up_dis": self.option("up_dis"),
            "down_dis": self.option("down_dis"),
        }
        self.target_predict_known.set_options(options)
        self.target_predict_known.run()

    def run_novol(self):
        options = {
            "lncrna_gtf" : self.option("novol"),
            "mrna_gtf" : self.option("mrna_gtf"),
            "annotation" : self.option("annotation"),
            "up_dis": self.option("up_dis"),
            "down_dis": self.option("down_dis"),

        }
        self.target_predict_novol.set_options(options)
        self.target_predict_novol.run()

    def run_merge_target(self):
        options = {
            "target_gtf" : self.option("mrna_gtf"),
            "known_lnc_gtf" : self.option("known"),
            "novel_lnc_gtf" : self.option("novol"),
            "cis_known" : self.target_predict_known.output_dir + "/lnc_rna_cistarget.annot.xls",
            "cis_novol" : self.target_predict_novol.output_dir + "/lnc_rna_cistarget.annot.xls",
            "annotation" : self.option("annotation"),
            "exp_corr": self.expcorr_tool.output_dir + "/corr.xls",

        }
        self.merge_target.set_options(options)
        self.merge_target.run()



    def run(self):
        super(TargetCistransModule, self).run()
        self.on_rely([self.expcorr_tool, self.target_predict_known, self.target_predict_novol], self.run_merge_target)
        self.merge_target.on('end', self.set_output)
        self.filter_exp()
        self.run_exp_corr()
        self.run_known()
        self.run_novol()


    def change_and_sort(self, file_in, file_out):
        a = pd.read_table(file_in, sep="\t", header=0, index_col=0)
        # a.lnc_type = a.lnc_type.replace({"known": "ref", "novel": "new"})
        # a.sort_values(by=['lnc_type'], ascending=[False], inplace=True)
        a.sort_values(by=['lnc_type'], ascending=[True], inplace=True)
        a.to_csv(file_out, sep="\t")

    def set_output(self):
        all_files = os.listdir(self.merge_target.output_dir)
        for each in all_files:
            if each.endswith('annot.xls'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                self.change_and_sort(os.path.join(self.merge_target.output_dir, each), link)

        self.set_db()


    def set_db(self):
        self.end()

    def end(self):

        target_dir = self.output_dir


        repaths = [
            [".", "", "靶基因"],
            ["known_lncrna_cistarget.annot.xls", "", "已知miRNA 对应的靶基因详情表"],
            ["novol_lncrna_cistarget.annot.xls", "", "新miRNA 对应的靶基因详情表"],
        ]

        super(TargetCistransModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''


    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet

        lnc_dir = "/mnt/ilustre/users/sanger-dev/workspace/20191129/WholeTranscriptome_tsg_36305/Transfer/remote_input/indir/workflow_results"

        known = os.path.join(lnc_dir,
                             'large_gush/filter_by_express/filtered_file/known_lncrna.gtf')
        novol = os.path.join(lnc_dir,
                             'large_gush/filter_by_express/filtered_file/novel_lncrna.gtf')
        mrna_gtf = os.path.join(lnc_dir,
                                'large_gush/filter_by_express/filtered_file/all_mrna.gtf')
        annotation = os.path.join(lnc_dir, 'annotation/allannot_class/all_annot.xls')
        up_dis = 20
        down_dis = 20
        exp_way = "tpm"
        exp_matrix_lnc = os.path.join(lnc_dir,
                                      'large_gush/filter_by_express/classifyquant/T.{}.txt'.format(exp_way))
        exp_matrix_target = os.path.join(lnc_dir,
                                         'large_gush/filter_by_express/classifyquant/G.{}.txt'.format(exp_way))
        opts = {
            'known': known,
            'novol': novol,
            'mrna_gtf': mrna_gtf,
            'annotation': annotation,
            'up_dis': up_dis,
            'down_dis': down_dis,
            'exp_matrix_lnc': exp_matrix_lnc,
            'exp_matrix_target': exp_matrix_target,
            'cor_cutoff': 0.8,
            'corr_way': "spearmanr",
            'padjust_way': "bh",
            'filter_lnc_file': '/mnt/ilustre/users/sanger-dev/workspace/20191129/WholeTranscriptome_tsg_36305/DiffSplit/output/lncrna/summary.txt',
            'filter_target_file': '/mnt/ilustre/users/sanger-dev/workspace/20191129/WholeTranscriptome_tsg_36305/Transfer/output/diff_exp_g/summary.txt'
        }


        data = {
            'id': 'target_cistrans_{}'.format(random.randint(1000, 9999)),
            'type': 'module',
            'name': 'whole_transcriptome.target_cistrans',
            'instant': False,
            'options': opts
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
