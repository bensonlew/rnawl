# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import json
import os
import shutil
import unittest

import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class CoxAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(CoxAgent, self).__init__(parent)
        options = [
            {'name': 'meta_file', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'exp_file', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'factor_list', 'type': 'string', },
            {'name': 'gene_list', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'cox_table', 'type': 'outfile', 'format': 'medical_transcriptome.common'},
            {'name': 'cox_graph', 'type': 'outfile', 'format': "medical_transcriptome.common"}

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
        super(CoxAgent, self).end()


class CoxTool(Tool):
    def __init__(self, config):
        super(CoxTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': 'bioinfo/miniconda2/bin/Rscript',
        }
        self.script = {
            'cox': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/cox.py'),
        }
        self.file = {
            'cox_table': os.path.join(self.output_dir, 'cox.txt'),
            'cox_graph': os.path.join(self.output_dir, 'cox.png')
        }

    def run(self):
        super(CoxTool, self).run()
        self.cox()
        self.set_output()
        self.end()

    def cox(self):
        cmd = '{} {}'.format(self.program['python'], self.script['cox'])
        cmd += ' -m {}'.format(self.option('meta_file').path)
        cmd += ' -e {}'.format(self.option('exp_file').path)
        cmd += ' -f {}'.format(self.option('factor_list'))
        cmd += ' -g {}'.format(self.option('gene_list').path)
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_cox'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33704402")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33704403")

    def set_output(self):
        self.option('cox_table').set_path(self.file['cox_table'])
        self.option('cox_graph').set_path(self.file['cox_graph'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'cox{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.tool.cox',
            'instant': False,
            'options': {
                'meta_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/lung_test.txt',
                'exp_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/count_test.txt',
                'factor_list': 'sex;age',
                'gene_list': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/gene'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
