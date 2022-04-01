# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class PrepareAgent(Agent):
    '''
    last_modify: 2019.04.10
    '''
    def __init__(self, parent):
        super(PrepareAgent, self).__init__(parent)
        options = [
            {'name': 'raw_fasta', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'txpt2gene', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 't2g', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'g2t', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'clean_fasta', 'type': 'outfile', 'format': 'ref_rna_v2.fasta'},
        ]
        self.add_option(options)
        self.step.add_steps('prepare')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.prepare.start()
        self.step.update()

    def stepfinish(self):
        self.step.prepare.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(os.path.getsize(self.option('raw_fasta').path)) / 1024 ** 3 + 4))

    def end(self):
        super(PrepareAgent, self).end()

class PrepareTool(Tool):
    def __init__(self, config):
        super(PrepareTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.prepare_py = os.path.join(self.config.PACKAGE_DIR, 'tool_lab/expression/prepare.py')
        self.t2g_pairs = os.path.join(self.output_dir, 't2g.pairs')
        self.g2t_pairs = os.path.join(self.output_dir, 'g2t.pairs')
        self.clean_fasta = os.path.join(self.output_dir, 'transcripts.fasta')

    def run(self):
        super(PrepareTool, self).run()
        self.run_prepare()
        self.set_output()
        self.end()

    def run_prepare(self):
        cmd = '{} {}'.format(self.python, self.prepare_py)
        cmd += ' -i {}'.format(self.option('raw_fasta').path)
        cmd += ' -m {}'.format(self.option('txpt2gene').path)
        cmd += ' -t {}'.format(self.t2g_pairs)
        cmd += ' -g {}'.format(self.g2t_pairs)
        cmd += ' -o {}'.format(self.clean_fasta)
        cmd_name = 'run_prepare'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('succeed in rerunning {}'.format(cmd_name))
            else:
                self.set_error('fail to rerun {}, abord'.format(cmd_name))
        else:
            self.set_error('fail to run {}, abord'.format(cmd_name))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('t2g').set_path(self.t2g_pairs)
        self.option('g2t').set_path(self.g2t_pairs)
        self.option('clean_fasta').set_path(self.clean_fasta)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'prepare_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.prepare',
            'instant': False,
            'options': {
                'raw_fasta': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'txpt2gene': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/trans2gene',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
