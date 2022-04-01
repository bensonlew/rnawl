# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from biocluster.workflow import Workflow
import os
import glob
import unittest
import json
import re
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class DiffexpBatchWorkflow(Workflow):
    '''
    差异分析
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffexpBatchWorkflow, self).__init__(wsheet_object)
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
            dict(name="filter_method", type="string", default=None),  # filter_method 20190801修改
            dict(name="tpm_filter_threshold", type="float", default="0"),  # 20190708 添加 by fwy
            dict(name="padjust_way", type='string', default="BH"),
            dict(name='exp_type', type='string', default='tpm'),
            # prob NOIseq独有的阈值
            dict(name='prob', type='float', default=0.8),
            # method: DESeq2, edgeR, DEGseq, NOIseq, Limma
            dict(name="method", type="string", default="DESeq2"),
            # batch by zjx 20200917
            dict(name="is_batch", type="bool", default=False),
            dict(name="has_batch", type="bool"),
            dict(name="batch_matrix", type="infile", format="denovo_rna_v2.common"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("small_rna_v2.diff.diffexp_batch")
        # self.all_exp = self.api.api("small_rna_v2.all_exp")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/02 Diff_Express')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(DiffexpBatchWorkflow, self).send_log(data)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for opt in self.sheet.options():
            self.logger.debug('{} - {}'.format(opt, self.option(opt)))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        '''
        define running logic
        '''
        self.tool.on("end", self.run_uniform)
        self.get_run_log()
        self.run_tool()
        super(DiffexpBatchWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_diff", main_id=self.option('diff_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        opts = dict(
            exp=self.option('exp_matrix'),
            exp_type=self.option('exp_type'),
            count=self.option('count'),
            group=self.option('group'),
            method=self.option('method'),
            cmp=self.option('cmp'),
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
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=self.option('count'),
                    group=self.option('group'),
                    method=self.option('method'),
                    cmp=self.option('cmp'),
                    pvalue=self.option('pvalue'),
                    pvalue_padjust=self.option('pvalue_padjust'),
                    padjust_way=self.option('padjust_way'),
                    fc=self.option('fc'),
                    tpm_filter_threshold=self.option("tpm_filter_threshold"),
                    is_batch=self.option('is_batch'),
                    has_batch=self.option('has_batch'),
                    batch_matrix=self.option('batch_matrix'),
                    filter_method=self.option("filter_method"),
                    no_filter=True,
                    analysis='split'
                )
            else:
                opts = dict(
                    exp=self.option('exp_matrix'),
                    exp_type=self.option('exp_type'),
                    count=self.option('count'),
                    group=self.option('group'),
                    method=self.option('method'),
                    cmp=self.option('cmp'),
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
        """
        保存结果表到mongo数据库中
        """
        self.stop_timeout_check()
        self.all_exp = self.api.api("small_rna_v2.all_exp")
        # add result info
        self.all_exp.add_diffexp_all(self.uniform.output_dir, self.tool.output_dir,
                                     group_dict=self.option('group_dict'),
                                     main_id=self.option('diff_main_id'),
                                     exp_id=self.option('exp_id'),
                                     diff_method=self.option('method'),
                                     create_geneset=False,
                                     pvalue_padjust=self.option('pvalue_padjust')
                                     )

        self.end()

    def end(self):
        '''
        declare output_dir and call super class method to end workflow
        '''
        # if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
        #     os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        # os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["02 Diff_Express", "", "miRNA表达量差异分析结果目录",0],
        ]
        result_dir.add_relpath_rules([
            ['.', '', 'miRNA表达量差异分析文件', 0],

        ])
        result_dir.add_regexp_rules([
            [r'.*summary.*\.xls', 'xls', '差异表达miRNA统计表', 0],
            [r'.*_vs_.*\.xls', '', '差异表达miRNA详情表（单个比较组）', 0],
            [r".*_vs_.*\..*normalize\.xls", "xls", "差异矩阵标准化结果表", 0],
            [r".*_vs_.*\..*sizeFactor\.xls", "xls", "差异矩阵标准化因子表", 0],
            [r".*_vs_.*\..*normFactor\.xls", "xls", "差异矩阵标准化因子表", 0],
            [r'.*total_diff_stat.*\.xls', 'xls', '差异表达miRNA详情表(所有比较组）', 0],
            [r'run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(DiffexpBatchWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):

        from mbio.workflows.small_rna_v2.report.diffexp import DiffexpWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna_v2.report.diffexp',
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/all_miR_norm.xls',
                'group_dict': json.dumps({'GH': ['GH1', 'GH2'], 'NFAs': ['NFAs1', 'NFAs2'],
                                          'Normal': ['Normal1', 'Normal2'], 'PRL': ['PRL1', 'PRL2']}),
                'count': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/all_miR_count.xls',
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
        wf.sheet.id = 'small_rna_v2'
        wf.sheet.project_sn = 'small_rna_v2'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)