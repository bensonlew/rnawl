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


class ExpPcaWorkflow(Workflow):
    '''
    TODO: description
    '''
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ExpPcaWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'exp_matrix', 'type': 'string'},
            {'name': 'group_dict', 'type': 'string'},  # to_file require, no use here
            {'name': 'pca_main_id', 'type': 'string'},
            {'name': 'ellipse', 'type': 'string'},
            {'name': 'group_table', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'update_info', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Express/02 Exp_PCA')
        self.inter_dirs = []
        self.tool = self.add_tool("small_rna.exp_pca")
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
        super(ExpPcaWorkflow, self).send_log(data)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for o in self.sheet.options():
            self.logger.debug('{} - {}'.format(o, self.option(o)))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def run(self):
        '''
        define running logic
        '''

        self.get_run_log()
        self.run_tool()
        super(ExpPcaWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_exp_pca", main_id=self.option('pca_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_tool(self):
        '''
        set options for tool
        set trigger point
        '''
        options = {
            'express_matrix': self.option('exp_matrix'),
        }
        if self.option('ellipse') == 'yes':
            self.tool.on("end", self.run_ellipse)
        else:
            self.tool.on("end", self.set_output)
        self.tool.set_options(options)
        self.tool.run()

    def run_ellipse(self):
        self.ellipse = self.add_tool('graph.ellipse')
        self.ellipse.set_options({
            'analysis': 'pca',
            'group_table': self.option('group_table').prop['path'],
            'pc_table': os.path.join(self.tool.output_dir, 'PCA.xls'),
        })
        self.ellipse.on('end', self.set_output)
        self.ellipse.run()

    def set_output(self):
        '''
        link result in output_dir of tool to output_dir of workflow
        '''
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for source in glob.glob(os.path.join(self.tool.output_dir, '*')):
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
        pca_output_dir = self.output_dir
        self.all_exp.add_exp_pca(pca_output_dir, main_id=self.option('pca_main_id'))
        if self.option('ellipse') == 'yes':
            self.all_exp.insert_ellipse_table(self.ellipse.work_dir + '/ellipse_out.xls',
                                           main_id=self.option('pca_main_id'))
        if self.option('ellipse') == 'no':
            pass
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
            ["01 Express/02 Exp_PCA", "", "样本间PCA分析", 0]
        ]
        result_dir.add_relpath_rules([
            ['.', '', '样本间PCA分析文件', 0],
            ['PCA.xls', '', '样本间PCA分析结果表', 0],
            ['Explained_variance_ratio.xls', '', '主成分解释表', 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(ExpPcaWorkflow, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test_add_known_tpm(self):

        from mbio.workflows.small_rna.report.exp_pca import ExpPcaWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            'id': 'workflow_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna.report.exp_pca',
            'options': {
                'exp_matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_norm.xls',
                'group_dict': json.dumps({'GH': ['GH1', 'GH2'], 'NFAs': ['NFAs1', 'NFAs2'],
                                          'Normal': ['Normal1', 'Normal2'], 'PRL': ['PRL1', 'PRL2']}),
                'pca_main_id': None
            }
        }

        wsheet = Sheet(data=data)
        wf = ExpPcaWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()