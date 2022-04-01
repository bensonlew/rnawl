# -*- coding: utf-8 -*-
# __author__ = "chenyanyan, 2016.10.12"
# last_modify by khl 20170504

from biocluster.workflow import Workflow
import os, re, glob
import pandas as pd
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import unittest


class CoxWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CoxWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'meta_file', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'exp_file', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_list', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'gene_str', 'type': 'string', 'default': None},
            {'name': 'sample_str', 'type': 'string', 'default': 'All'},
            {'name': 'is_show', 'type': 'string'}, #没有用的参数
            {'name': 'gene_choose', 'type': 'string'}, #没有用的参数
            {'name': 'cox_table', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'cox_graph', 'type': 'outfile', 'format': "ref_rna_v2.common"},
            {'name': 'update_info', 'type': 'string'},
            {'name': "main_id", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.cox = self.add_tool("tool_lab.cox")

    def run(self):
        self.run_cox()
        super(CoxWorkflow, self).run()

    def run_cox(self):
        with open(self.option('meta_file').path, 'r') as m:
            a = m.readline()
            factor_list = a.strip().split('\t')[3:]
            factor_str = ';'.join(factor_list)
        meta_file_filter_path = os.path.join(self.work_dir, 'meta_file.txt')
        exp_file_filter_path = os.path.join(self.work_dir, 'exp_file.txt')
        if self.option('sample_str').lower() != 'all':
            meta_file = pd.read_table(self.option('meta_file').path, header=0, index_col=0, sep='\t')
            meta_file.index = meta_file.index.astype('str', copy=False)  # added by zhangyitong on 20211011
            sample_list = self.option('sample_str').split(';')
            meta_file_filter = meta_file.loc[sample_list]
            meta_file_filter.to_csv(meta_file_filter_path, header=True, index=True, sep='\t')
            if self.option('exp_file').is_set:
                exp_file = pd.read_table(self.option('exp_file').path, header=0, index_col=0, sep='\t')
                exp_file_filter = exp_file[sample_list]
                exp_file_filter.fillna(0, inplace=True)     # added by zhangyitong on 20211012
                if (exp_file_filter < 0).values.any():
                    self.set_error('输入的表达谱存在负值，请检查。')
                exp_file_filter.to_csv(exp_file_filter_path, header=True, index=True, sep='\t')
            else:
                exp_file_filter_path = None
        else:
            if self.option('exp_file').is_set:
                exp_file = pd.read_table(self.option('exp_file').path, header=0, index_col=0, sep='\t')
                sample_exp = exp_file.columns.tolist()
                meta_file = pd.read_table(self.option('meta_file').path, header=0, index_col=0, sep='\t')
                meta_file.index = meta_file.index.astype('str', copy=False)
                sample_meta = meta_file.index.tolist()
                intersection_list = list(set(sample_exp).intersection(set(sample_meta)))
                meta_file_filter = meta_file.loc[intersection_list]
                meta_file_filter.to_csv(meta_file_filter_path, header=True, index=True, sep='\t')
                exp_file_filter = exp_file[intersection_list]
                exp_file_filter.fillna(0, inplace=True)
                if (exp_file_filter < 0).values.any():
                    self.set_error('输入的表达谱存在负值，请检查。')
                exp_file_filter.to_csv(exp_file_filter_path, header=True, index=True, sep='\t')
            else:
                meta_file = pd.read_table(self.option('meta_file').path, header=0, index_col=0, sep='\t')
                meta_file.to_csv(meta_file_filter_path, header=True, index=True, sep='\t')
        if self.option('gene_str'):
            gene_list_path = os.path.join(self.work_dir, 'gene_list')
            with open(gene_list_path, 'w') as g:
                g.write('seq_id' + '\n')
                gene_list = self.option('gene_str').split(';')
                for i in gene_list:
                    g.write(i + '\n')
            self.option('gene_list').set_path(gene_list_path)
        if self.option('exp_file').is_set:
            opts = {
                'meta_file': meta_file_filter_path,
                'exp_file': exp_file_filter_path,
                'factor_list': factor_str,
                'gene_list': self.option('gene_list'),
            }
        else:
            opts = {
                'meta_file': meta_file_filter_path,
                'factor_list': factor_str,
            }
        self.cox.set_options(opts)
        self.cox.on('end', self.set_db)
        self.cox.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        pdf_path = self.cox.option('cox_graph').path
        s3_output = '{}/{}'.format(self._sheet.output, os.path.basename(pdf_path))
        tool_cox = self.api.api("tool_lab.tool_cox")
        # # add result info
        tool_cox.add_cox(self.option('main_id'), self.cox.option('cox_table').path, s3_output)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.cox.output_dir)
        # self.inter_dirs = [
        #     ["04 GeneSet", "", "基因集分析结果目录",0],
        #     ["04 GeneSet/01 Cluster_Analysis", "", "聚类分析", 0]
        # ]
        # result_dir.add_relpath_rules([
        #     [".", "", "聚类分析文件",0,"211530"],
        #     ["./seq.cluster_tree.txt", "txt", "基因/转录本聚类树文件",0,"211531"],
        #     ["./seq.kmeans_cluster.txt", "txt", "基因/转录本聚类树文件", 0],
        #     ["./sample.cluster_tree.txt", "txt", "样本聚类树文件",0,"211532"],
        #     ["./expression_matrix.xls", "xls", "聚类热图分析表",0,"211533"],
        #     ["./heatmap.pdf", 'pdf', "聚类热图",0],
        # ])
        # result_dir.add_regexp_rules([
        #     [r'.*subcluster_.*\.xls', 'xls', '子聚类分析表',0,"211534"],
        # ])
        super(CoxWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.medical_transcriptome.report.cox import CoxWorkflow
        from biocluster.wsheet import Sheet
        import random
        # test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "estimate" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "medical_transcriptome.cox",
            "output":"s3://medical_transcriptome/files/test/medical_transcriptome/medical_transcriptome/interaction_results/ASprofile_20200813_092731",
            "options": {
                'exp_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/count_test.txt',
            }
        }

        wsheet = Sheet(data=data)
        wf =CoxWorkflow(wsheet)
        wf.sheet.id = 'medical_transcriptome'
        wf.sheet.project_sn = 'medical_transcriptome'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)