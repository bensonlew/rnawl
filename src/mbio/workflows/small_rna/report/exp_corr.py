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


class ExpCorrWorkflow(Workflow):
    '''
    transfer exp_matrix (file:str), scm (str), scd (str) and corr_method (str) to tool
    obtain output_dir from tool
    deliver output_dir to api according to main_id of sg_exp_venn
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpCorrWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'string'},
            {'name': 'group_dict', 'type': 'string'}, # to_file require, no use here
            {'name': 'use_log', 'type': 'bool', 'default': False}, # to_file require, no use here
            {'name': 'corr_method', 'type': 'string', 'default': 'pearson'},
            {'name': 'scm', 'type': 'string', 'default': 'complete'},
            {'name': 'scd', 'type': 'string', 'default': 'euclidean'},
            {'name': 'corr_main_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Express/01 Exp_Corr')
        self.inter_dirs = []
        self.tool = self.add_tool('small_rna.exp_corr')
        self.all_exp = self.api.api('small_rna.all_exp')

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
        super(ExpCorrWorkflow, self).send_log(data)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for o in self.sheet.options():
            self.logger.debug('{} - {}'.format(o, self.option(o)))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        '''
        define running logic
        '''
        self.tool.on('end', self.set_output)
        self.get_run_log()
        self.run_tool()
        super(ExpCorrWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_exp_corr", main_id=self.option('corr_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        '''
        set options for tool
        set trigger point
        '''
        options = {
            'express_matrix': self.option('exp_matrix'),
            'scm': self.option('scm'),
            'scd': self.option('scd'),
            'corr_method': self.option('corr_method'),
        }
        self.tool.set_options(options)
        self.tool.run()

    def set_output(self):
        '''
        link result in output_dir of tool to output_dir of workflow
        '''
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = os.path.join(self.tool.output_dir, 'sample_correlation.xls')
        link_name = os.path.join(self.output_dir, 'sample_correlation.xls')
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.set_db()

    def set_db(self):
        '''
        link result in output_dir of tool to output_dir of workflow
        '''
        self.logger.info('start set_db at {}'.format(self.__class__.__name__))
        corr_output_dir = self.tool.output_dir
        self.all_exp.add_exp_corr(corr_output_dir, main_id=self.option('corr_main_id'))
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
        self.inter_dirs = [
            ["01 Express", "", "miRNA表达量分析结果目录",0],
            ["01 Express/01 Exp_Corr", "", "样本间相关性分析", 0]
        ]
        result_dir.add_relpath_rules([
            ['.', '', '样本间相关性分析文件', 0],
            ['sample_correlation.xls', '', '样本间相关性系数表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(ExpCorrWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_add_known_tpm(self):

        from mbio.workflows.small_rna.report.exp_corr import ExpCorrWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna.report.exp_corr',
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_norm.xls',
                'group_dict': json.dumps({'GH': ['GH1', 'GH2'], 'NFAs': ['NFAs1', 'NFAs2'],
                                          'Normal': ['Normal1', 'Normal2'], 'PRL': ['PRL1', 'PRL2']}),
                'scm': 'complete',
                'scd': 'correlation',
                'corr_method': 'pearson',
                'corr_main_id': None
            }
        }

        wsheet = Sheet(data=data)
        wf = ExpCorrWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()