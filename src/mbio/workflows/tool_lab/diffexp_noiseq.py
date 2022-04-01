# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import pandas as pd
import glob
import os
import json
import time
import re
from biocluster.file import getsize, exists
from biocluster.file import download
import pandas as pd
import unittest


class DiffexpNoiseqWorkflow(Workflow):
    """
    差异分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpNoiseqWorkflow, self).__init__(wsheet_object)
        options = [

            dict(name="main_id", type="string"),
            dict(name="count", type="infile", format="ref_rna_v2.common"),
            dict(name="group", type="infile", format='ref_rna_v2.common'),
            dict(name="cmp", type="infile", format='ref_rna_v2.common'),
            dict(name="task_id", type="string"),
            dict(name="is_batch", type="bool", default=False),  #是否做批次效应处理
            dict(name="has_batch", type="bool", default=True),  #是否有批次处理表
            # dict(name="batch_matrix", type="infile", format="ref_rna_v2.common"),  #上传批次效应表
            # pvalue 统计判断阈值
            # dict(name="pvalue", type="float", default=0.05),
            # dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            # dict(name="padjust_way", type='string', default="BH"),
            # method: DESeq2, edgeR, DEGseq, limma, NOIseq
            dict(name="method", type="string"),
            # dict(name='prob', type='float', default=0.8),
            dict(name='norm', type='string', default='TMM'),
            # dict(name='model', type='string', default='glmQLFit'),
            # dict(name='pair_table', type='infile', format='ref_rna_v2.common'),
            dict(name='prob',type='float', default=0.9),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("tool_lab.diffexp")
        self.diffexp = self.api.api("tool_lab.diffexp")

    def prepare(self):
        count = pd.read_table(self.option('count').path, header=0, index_col=0, sep='\t')
        count.index.name = 'seq_id'
        count.to_csv(self.option('count').path, header=True, index=True, sep='\t')

    def run(self):
        self.tool.on("end", self.set_db)
        self.prepare()
        self.run_tool()
        super(DiffexpNoiseqWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        output_dir = self.tool.output_dir
        pdf_path = glob.glob(os.path.join(output_dir, '*.pdf'))
        pdf_list = list()
        pdf_name = list()
        for i in pdf_path:
            pdf_list.append('{}/{}'.format(self._sheet.output, os.path.basename(i)))
            pdf_name.append(os.path.basename(i))
        # add result info

        self.diffexp.add_diffexp_noiseq(self.tool.output_dir,
            main_id=self.option('main_id'),
            diff_method='noiseq',
            create_geneset=False,
            s3_output=pdf_list,
            pdf_name=pdf_name
            )
        # self.paste_annotation()
        self.set_output()
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "差异分析结果目录", 0, "211063"],
        ])
        result_dir.add_regexp_rules([
            [r'.*_vs_.*\.xls', 'xls', '差异表达基因详情表',0,"211518"],
            [r'.*summary.*\.xls', 'xls', '差异表达基因统计表',0,"211519"],
            [r'.*total_diff_stat.*\.xls', 'xls', '差异表达基因详情总表',0,"211520"],
        ])
        super(DiffexpNoiseqWorkflow, self).end()

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
    def run_tool(self):
        opts = dict(
            #基本参数
            count=self.option('count').path,
            group=self.option('group').path,
            method='NOIseq',
            cmp=self.option('cmp').path,
            fc=self.option('fc'),
            prob=self.option('prob'),
            # pvalue=self.option('pvalue'),
            # pvalue_padjust=self.option('pvalue_padjust'),
            # padjust_way=self.option('padjust_way'),
            #高级参数
            norm=self.option('norm'),
            # model=self.option('model')


        )
        # if self.option('pair_table').is_set:
        #     opts.update({'is_batch': True, 'has_batch': True, 'batch_matrix': self.option('pair_table').path})




        self.tool.set_options(opts)
        self.tool.run()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.diffexp_noiseq import DiffexpNoiseqWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            'id': 'diff_exp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'tool_lab.diffexp_noiseq',
            "options": dict(
                count='/mnt/ilustre/users/sanger-dev/workspace/20200629/Refrna_tsg_37857/Quant/output/gene.count.matrix',
                group="/mnt/ilustre/users/sanger-dev/workspace/20200629/Refrna_tsg_37857/remote_input/group_table/group.txt",
                cmp= "/mnt/ilustre/users/sanger-dev/workspace/20200629/Refrna_tsg_37857/remote_input/control_file/compare.txt",
                # padjust_way="BH",
                # pvalue=0.05,
                fc=2,
                prob=0.9,
                norm='TMM'


            )
        }
        wsheet = Sheet(data=data)
        wf = DiffexpNoiseqWorkflow(wsheet)
        wf.sheet.project_sn = '188_5dba6f542345b'
        wf.sheet.task_id = '123456'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


