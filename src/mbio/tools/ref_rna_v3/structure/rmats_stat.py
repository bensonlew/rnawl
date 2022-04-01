# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.ref_rna_v3.functions import toolfuncdeco, runcmd


class RmatsStatAgent(Agent):
    '''
    last_modify: 2019.07.17
    '''

    def __init__(self, parent):
        super(RmatsStatAgent, self).__init__(parent)
        options = [
            {'name': 'root', 'type': 'infile', 'format': 'ref_rna_v3.common_dir'},
            {'name': 'pvalue_fdr', 'type': 'string', 'default': None},
            {'name': 'fdr', 'type': 'float', 'default': None},
            {'name': 'psi', 'type': 'float', 'default': None},
            {'name': 'main_id', 'type': 'string', 'default': None},
            {'name': 'event_stats', 'type': 'outfile', 'format': 'ref_rna_v3.common'},
            {'name': 'psi_stats', 'type': 'outfile', 'format': 'ref_rna_v3.common'},
        ]
        self.add_option(options)
        self.step.add_steps('rmats_stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def step_start(self):
        self.step.rmats_stat.start()
        self.step.update()

    def step_end(self):
        self.step.rmats_stat.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    @toolfuncdeco
    def end(self):
        super(RmatsStatAgent, self).end()


class RmatsStatTool(Tool):
    def __init__(self, config):
        super(RmatsStatTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'rmats_stats': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/structure/rmats_stat.py')
        }
        self.file = {
            'event_stats': os.path.join(self.output_dir, 'event_stats.file.txt'),
            'psi_stats': os.path.join(self.output_dir, 'psi_stats.file.txt')
        }

    @toolfuncdeco
    def run(self):
        super(RmatsStatTool, self).run()
        self.run_rmats_stat()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_rmats_stat(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats_stats'])
        cmd += ' -i {}'.format(self.option('root').path)
        cmd += ' -m {}'.format(self.option('pvalue_fdr'))
        cmd += ' -s {}'.format(self.option('main_id'))
        cmd += ' -c {}'.format(self.option('fdr'))
        cmd += ' -p {}'.format(self.option('psi'))
        cmd += ' -o {}'.format(self.output_dir)
        if self.config.DBVersion == 1:
            DBVersion =1
        else:
            DBVersion = 0
        cmd += ' -v {}'.format(DBVersion)
        cmd_name = 'run_rmats_stat'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def set_output(self):
        self.option('event_stats').set_path(self.file['event_stats'])
        self.option('psi_stats').set_path(self.file['psi_stats'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_stat_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.structure.rmats_stat',
            'instant': False,
            'options': {
                'root': '/mnt/ilustre/users/sanger-dev/workspace/20190617/RmatsStat_tsg_33555_4445_6310/remote_input/root/Con_34_vs_Vit_34',
                'pvalue_fdr': 'fdr',
                'fdr': 0.05,
                'psi': 0.1
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
