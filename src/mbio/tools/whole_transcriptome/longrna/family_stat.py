# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class FamilyStatAgent(Agent):
    '''
    last_modify: 2019.04.19
    '''
    def __init__(self, parent):
        super(FamilyStatAgent, self).__init__(parent)
        options = [
            {'name': 'tsv_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'known_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'novel_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 't2g', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'tabular', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('family_stat')
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.family_stat.start()
        self.step.update()

    def step_finish(self):
        self.step.family_stat.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(FamilyStatAgent, self).end()

class FamilyStatTool(Tool):
    def __init__(self, config):
        super(FamilyStatTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.family_statistics_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/family_statistics.py')
        self.tabular = os.path.join(self.output_dir, 'lncRNA_family.tabular')

    def run(self):
        super(FamilyStatTool, self).run()
        self.run_family_stat()
        self.set_output()
        self.end()

    def run_family_stat(self):
        cmd = '{} {}'.format(self.python, self.family_statistics_py)
        cmd += ' -i {}'.format(self.option('tsv_list').path)
        cmd += ' -k {}'.format(self.option('known_list').path)
        cmd += ' -n {}'.format(self.option('novel_list').path)
        cmd += ' -t {}'.format(self.option('t2g').path)
        cmd += ' -o {}'.format(self.tabular)
        cmd_name = 'run_family_stat'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False, block=True):
        if shell:
            cmd = self.config.SOFTWARE_DIR + '/' + cmd
        command = self.add_command(cmd_name, cmd, shell=shell)
        command.run()
        command.no_check = True
        if block:
            self.wait()
            for n, c in self.commands.items():
                if c.no_check:
                    if c.return_code == c.default_return_code:
                        c.no_check = False
                        self.logger.info('succeed in running {}'.format(n))
                    else:
                        self.set_error('fail to run {}, abord'.format(n))

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        self.option('tabular').set_path(self.tabular)
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
            'id': 'family_stat_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.family_stat',
            'instant': False,
            'options': {
                'tsv_list': '/mnt/ilustre/users/sanger-dev/workspace/20190419/Single_lncrna_family_1515_6581/LncrnaFamily/tsv.list',
                'known_list': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33905/FilterByExpress/output/filtered_file/known_lncrna_ids.list',
                'novel_list': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33905/FilterByExpress/output/filtered_file/novel_lncrna_ids.list',
                't2g': '/mnt/ilustre/users/sanger-dev/workspace/20190418/LncRna_tsg_33905/Assemble/output/NewTranscripts/trans2gene',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
