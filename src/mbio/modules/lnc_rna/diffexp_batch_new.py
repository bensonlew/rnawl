# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import pandas as pd
import os
import glob
from biocluster.file import getsize, exists
from biocluster.file import download
from mbio.packages.ref_rna_v2.functions import tryforgood
import unittest
from biocluster.module import Module


class DiffexpBatchNewModule(Module):
    """
    差异分析
    """
    def __init__(self, work_id):
        super(DiffexpBatchNewModule, self).__init__(work_id)
        options = [
            dict(name='exp', type='infile', format='lnc_rna.express_matrix'),
            dict(name="group_dict", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name="count", type="infile", format='lnc_rna.express_matrix'),
            dict(name="group", type="infile", format='sample.group_table'),
            dict(name="cmp", type="infile", format='sample.control_table'),
            dict(name='exp_type', type='string', default="tpm"),
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="filter_method", type="string", default=None),  # filter_method 20190801修改
            dict(name="tpm_filter_threshold", type="float", default="0"), #20190708 添加 by fwy
            dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),
            # DESeq2 DE test method, Wald|LRT
            dict(name="deseq2_method", type="string", default="Wald"),
            # edger_method DE test method,exactTest|glmLRT|glmQLFTest
            dict(name="edger_method", type="string", default="glmQLFTest"),
            # degseq_method DE test method, LRT|CTR|FET|MARS|MATR|FC
            dict(name="degseq_method", type="string", default="MARS"),
            dict(name="task_id", type="string", default=""),
            dict(name="exp_level", type="string", default=""),
            # batch by zjx 20200917
            dict(name="is_batch", type="bool", default=False),
            dict(name="has_batch", type="bool"),
            dict(name="batch_matrix", type="infile", format="lnc_rna.common"),
            # noiseq
            dict(name='prob', type='float', default=0.8),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.tool = self.add_tool("lnc_rna.diffexp_batch")

    def run(self):
        super(DiffexpBatchNewModule, self).run()
        self.run_tool()

    def set_output(self):
        for file_name in os.listdir(self.tool.output_dir):
            source = os.path.join(self.tool.output_dir, file_name)
            link_name = os.path.join(self.output_dir, file_name)
            if os.path.isfile(link_name):
                os.remove(link_name)
            os.link(source, link_name)
        self.end()

    def end(self):
        super(DiffexpBatchNewModule, self).end()

    def run_tool(self):
        opts = dict(
            exp=self.option('exp').path,
            exp_type=self.option('exp_type'),
            count=self.option('count').path,
            group=self.option('group').path,
            method=self.option('method'),
            cmp=self.option('cmp').path,
            fc=self.option('fc'),
            tpm_filter_threshold=self.option("tpm_filter_threshold"),
            is_batch=self.option('is_batch'),
            filter_method=self.option("filter_method"),
            no_filter=True,
            analysis='split'
        )

        if self.option('is_batch') == False:
            if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
                opts.update(
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                )
            else:
                opts.update(
                    prob=self.option('prob')
                )
        else:
            if self.option('has_batch') == True:
                opts = dict(
                    exp=self.option('exp').path,
                    exp_type=self.option('exp_type'),
                    count=self.option('count').path,
                    group=self.option('group').path,
                    method=self.option('method'),
                    cmp=self.option('cmp').path,
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                    fc=self.option('fc'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    is_batch=self.option('is_batch'),
                    has_batch=self.option('has_batch'),
                    batch_matrix=self.option('batch_matrix').path,
                    filter_method=self.option("filter_method"),
                    no_filter=True,
                    analysis='split'
                )
            else:
                opts = dict(
                    exp=self.option('exp').path,
                    exp_type=self.option('exp_type'),
                    count=self.option('count').path,
                    group=self.option('group').path,
                    method=self.option('method').path,
                    cmp=self.option('cmp').path,
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                    fc=self.option('fc'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    is_batch=self.option('is_batch'),
                    has_batch=self.option('has_batch'),
                    filter_method=self.option("filter_method"),
                    no_filter=True,
                    analysis='split'
                )
        self.tool.set_options(opts)
        self.tool.on('end', self.run_uniform)
        self.tool.run()

    def run_uniform(self):
        self.uniform = self.add_tool('denovo_rna_v3.batch.uniform')
        opts = {
            'method': self.option('method'),
            'input_dir': self.tool.output_dir,
        }
        if self.option('method').lower() in ['deseq2', 'degseq', 'edger', 'limma']:
            opts.update({'pvalue_padjust': self.option('pvalue_padjust')})
        self.uniform.set_options(opts)
        self.uniform.on('end', self.set_output)
        self.uniform.run()



class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.denovo_rna_v3.report.diffexp_batch import DiffexpBatchWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'Diff_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'denovo_rna_v3.report.diffexp_batch',
            'options': dict(
                count="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/unigene.count.matrix.xls",
                exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/exp_matrix",
                method="DESeq2",
                group="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/group",
                cmp="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/cmp_little",
                padjust_way="BH",
                pvalue=0.05,
                fc=1,
                is_batch=False,
                exp_level='G',


            )
        }
        wsheet = Sheet(data=data)
        wf =DiffexpBatchWorkflow(wsheet)
        wf.sheet.id = 'diff'
        wf.sheet.project_sn = 'diff'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


