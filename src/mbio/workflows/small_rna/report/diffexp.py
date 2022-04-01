# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from biocluster.workflow import Workflow
import os
import glob
import unittest
import json
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class DiffexpWorkflow(Workflow):
    '''
    TODO: description
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='exp_matrix', type='string'),
            dict(name="group_dict", type="string"),
            dict(name="seq_type", type="string"),
            dict(name="diff_main_id", type="string"),
            dict(name="count", type="string"),
            dict(name="exp_id", type="string"),
            dict(name="group", type="string"),
            dict(name="cmp", type="string"),
            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="padjust_way", type='string', default="BH"),
            # prob NOIseq独有的阈值
            dict(name='prob', type='float', default=0.8),
            # method: DESeq2, edgeR, DEGseq, NOIseq, Limma
            dict(name="method", type="string", default="DESeq2"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("small_rna.diffexp")
        self.all_exp = self.api.api("small_rna.all_exp")

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for opt in self.sheet.options():
            self.logger.debug('{} - {}'.format(opt, self.option(opt)))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        '''
        define running logic
        '''
        self.tool.on("end", self.set_output)
        self.get_run_log()
        self.run_tool()
        super(DiffexpWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_diff", main_id=self.option('diff_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        '''
        set options for tool
        set trigger point
        '''
        options = dict(
            exp=self.option('exp_matrix'),
            count=self.option('count'),
            group=self.option('group'),
            method=self.option('method'),
            cmp=self.option('cmp'),
            # pvalue=self.option('pvalue'),
            # pvalue_padjust=self.option('pvalue_padjust'),
            # padjust_way=self.option('padjust_way'),
            fc=self.option('fc'),
        )
        if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma']:
            options.update(
                pvalue=self.option('pvalue'),
                pvalue_padjust=self.option('pvalue_padjust'),
                padjust_way=self.option('padjust_way'),
            )
        if self.option('method').lower() in ['noiseq']:
            options.update(
                prob=self.option('prob')
            )
        self.tool.set_options(options)
        self.tool.run()

    def set_output(self):
        '''
        link result in output_dir of tool to output_dir of workflow
        '''
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for source in glob.glob(os.path.join(self.tool.output_dir, '*.xls')):
            basename = os.path.basename(source)
            link_name = os.path.join(self.output_dir, basename)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        '''
        link result in output_dir of tool to output_dir of workflow
        '''
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        diff_output_dir = self.output_dir
        if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma']:
            self.all_exp.add_diffexp(diff_output=diff_output_dir,
                                     main_id=self.option('diff_main_id'),
                                     exp_id=self.option('exp_id'),
                                     diff_method=self.option('method'),
                                     create_geneset=False,
                                     pvalue_padjust=self.option('pvalue_padjust'))
        if self.option('method').lower() in ['noiseq']:
            self.all_exp.add_diffexp_noiseq(self.tool.output_dir,
                main_id=self.option('diff_main_id'),
                diff_method=self.option('method'),
                create_geneset=False,
            )
        self.logger.info('finish set_db at {}'.format(self.__class__.__name__))
        self.end()

    def end(self):
        '''
        declare output_dir and call super class method to end workflow
        '''
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', '表达量差异分析结果文件', 0],
            ['DESeq2_diff_summary.xls', '', '表达量差异统计表（所有比较组）', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'.*_vs_.*\.xls', '', '差异分析详情表（单个比较组）', 0],
            [r".*_vs_.*\..*normalize\.xls", "xls", "差异矩阵标准化结果表", 0],
            [r".*_vs_.*\..*sizeFactor\.xls", "xls", "差异矩阵标准化因子表", 0],
            [r".*_vs_.*\..*normFactor\.xls", "xls", "差异矩阵标准化因子表", 0],
        ])
        super(DiffexpWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_add_known(self):

        from mbio.workflows.small_rna.report.diffexp import DiffexpWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna.report.diffexp',
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_norm.xls',
                'group_dict': json.dumps({'GH': ['GH1', 'GH2'], 'NFAs': ['NFAs1', 'NFAs2'],
                                          'Normal': ['Normal1', 'Normal2'], 'PRL': ['PRL1', 'PRL2']}),
                'count': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_count.xls',
                'group': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/group.txt',
                'cmp': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/control.txt',
                'diff_main_id': None,
                'method': 'DESeq2',
                'fc': 2.0,
                'pvalue': 0.05,
                'pvalue_padjust': 'padjust',
                'padjust_way': 'BH',
            }
        }

        wsheet = Sheet(data=data)
        wf = DiffexpWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()