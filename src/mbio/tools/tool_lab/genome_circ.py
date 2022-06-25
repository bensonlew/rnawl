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


class GenomeCircAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(GenomeCircAgent, self).__init__(parent)
        options = [
            {'name': 'genome_type', 'type': 'string', 'default': 'finish'},
            {'name': 'genome_fa', 'type': 'string'},
            {'name': 'genome_struction', 'type': 'string'},
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
        super(GenomeCircAgent, self).end()


class GenomeCircTool(Tool):
    def __init__(self, config):
        super(GenomeCircTool, self).__init__(config)
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
            'genome_circ': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/genome_circ.py'),
        }
        self.file = {
            'output_file': os.path.join(self.output_dir, 'genome_circ.txt'),
        }

    def run(self):
        super(GenomeCircTool, self).run()
        self.run_genome_circ()
        self.set_output()
        self.end()

    def run_genome_circ(self):
        cmd = '{} {}'.format(self.program['python'], self.script['genome_circ'])
        cmd += ' -f {}'.format(self.option('genome_fa'))
        cmd += ' -s {}'.format(self.option('genome_struction'))
        cmd += ' -t {}'.format(self.option('genome_type'))
        cmd += ' -o {}'.format(self.file['output_file'])
        cmd_name = 'run_genome_circ'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704403")

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
