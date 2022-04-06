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


class PacbioAsAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(PacbioAsAgent, self).__init__(parent)
        options = [
            {'name': 'flnc', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        super(PacbioAsAgent, self).end()


class PacbioAsTool(Tool):
    def __init__(self, config):
        super(PacbioAsTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.python = os.path.join(software_dir, 'program/Python/bin')
        self.r = os.path.join(software_dir, 'program/R-3.3.1/bin')
        self.set_environ(PATH=self.python)
        self.set_environ(PATH=self.r)
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': os.path.join(software_dir, 'program/R-3.3.1/bin/Rscript'),
            'makeblastdb': "bioinfo/align/ncbi-blast-2.3.0+/bin/makeblastdb",
            'blastn': "bioinfo/align/ncbi-blast-2.3.0+/bin/blastn"
        }
        self.script = {
            'as': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/pacbio_as.py'),
        }
        self.file = {
            'ASdb': os.path.join(self.work_dir, 'ASdb'),
            'AS_blast': os.path.join(self.work_dir, 'AS_blast.out'),
            'output_file': os.path.join(self.output_dir, 'go_circ.pdf'),
        }

    def run(self):
        super(PacbioAsTool, self).run()
        self.run_as()
        self.run_blast()
        self.run_stat()
        self.set_output()
        self.end()

    def run_as(self):
        cmd = '{}'.format(self.program['makeblastdb'])
        cmd += ' -in {}'.format(self.option('flnc'))
        cmd += ' -dbtype nucl'
        cmd += ' -out ASdb -parse_seqids'
        cmd_name = 'run_as'
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

    def run_blast(self):
        'blastn -query {input_fasta} -out AS_blast.out -db ASdb -outfmt 5 -evalue 1e-5 -num_threads 4'
        cmd = '{}'.format(self.program['blastn'])
        cmd += ' -query {}'.format(self.option('flnc'))
        cmd += ' -out AS_blast.out'
        cmd += ' -db ASdb'
        cmd += ' -outfmt 5'
        cmd += ' -evalue 1e-5'
        cmd += ' -num_threads 4'
        cmd_name = 'run_blast'
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

    def run_stat(self):
        cmd = '{} {}'.format(self.program['python'], self.script['as'])
        cmd += ' -i {}'.format(self.file['AS_blast'])
        cmd += ' -o {}'.format(self.work_dir)
        cmd_name = 'run_stat'
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
