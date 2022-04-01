# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class SplitBedAgent(Agent):
    '''
    last_modify: 2019.02.27
    '''
    def __init__(self, parent):
        super(SplitBedAgent, self).__init__(parent)
        options = [
            {'name': 'ref_fa', 'type': 'infile', 'format':'lnc_rna.fasta'},
        ]
        self.add_option(options)
        self.step.add_steps('split_bed')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.split_bed.start()
        self.step.update()

    def step_end(self):
        self.step.split_bed.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('ref_fa').is_set:
            self.logger.debug('{} = {}'.format('ref_fa', self.option('ref_fa').prop['path']))
            self.infile_size = os.path.getsize(self.option('ref_fa').prop['path'])
        else:
            raise OptionError('reference FASTA file must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 4 + 4))

    def end(self):
        super(SplitBedAgent, self).end()

class SplitBedTool(Tool):
    def __init__(self, config):
        super(SplitBedTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.split_bed_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/split_bed.py')

    def run(self):
        super(SplitBedTool, self).run()
        self.run_split_bed()
        self.set_output()
        self.end()

    def run_split_bed(self):
        cmd = '{} {}'.format(self.python, self.split_bed_py)
        cmd += ' -i {}'.format(self.option('ref_fa').prop['path'])
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_split_bed'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
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
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_hsa(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'split_bed_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.split_bed',
            'instant': False,
            'options': {
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()