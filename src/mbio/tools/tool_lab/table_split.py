# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import json
import os
import shutil
import unittest
import math
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class TableSplitAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(TableSplitAgent, self).__init__(parent)
        options = [
            {'name': 'excel_file', 'type': 'string'},
            {'name': 'split_num', 'type': 'int'}
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(TableSplitAgent, self).end()


class TableSplitTool(Tool):
    def __init__(self, config):
        super(TableSplitTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.python = os.path.join(software_dir, 'miniconda2/bin')
        self.r = os.path.join(software_dir, 'program/R-3.3.1/bin')
        self.set_environ(PATH=self.python)
        self.set_environ(PATH=self.r)
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': os.path.join(software_dir, 'program/R-3.3.1/bin/Rscript'),
        }
        self.script = {
            'go_circ': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/GO_circ.py'),
        }
        self.file = {
            'output_file': os.path.join(self.output_dir, 'go_circ.pdf'),
        }

    def run(self):
        super(TableSplitTool, self).run()
        self.run_excel_split()
        self.set_output()
        self.end()

    def run_excel_split(self):
        file = pd.read_excel(self.option('excel_file'), header=0, sep='\t')
        columns_list = file.columns.tolist()
        split_name = columns_list[self.option('split_num')-1]
        file_list = file.groupby(split_name)
        for i in file_list:
            new_name = str(i[0]) + '_' +os.path.basename(self.option('excel_file'))
            file_name = os.path.join(self.output_dir, new_name)
            df_file = i[1]
            output_file = os.path.join(self.output_dir, file_name)
            df_file.to_csv(output_file, header=True, index=False, sep='\t')

    def set_output(self):
        # self.option('estimate_score_gct').set_path(self.file['estimate_score_gct'])
        # self.option('cluster_matrix').set_path(self.file['cluster_matrix'])
        pass
class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'estimate{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.tool.estimate',
            'instant': False,
            'options': {
                'exp_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/count_test.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
