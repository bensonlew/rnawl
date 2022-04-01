# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
import os
from biocluster.tool import Tool
import unittest

class Bowtie2BuildAgent(Agent):
    def __init__(self, parent):
        super(Bowtie2BuildAgent, self).__init__(parent)
        options = [
            {'name': 'reference_in', 'type': 'infile', 'format': 'ref_genome_db.fasta'},
            {'name': 'bt2_index_base', 'type': 'string', 'default': ''},
            {'name': 'threads', 'type': 'int', 'default': 8},
        ]
        self.add_option(options)
        self._memory_increase_step = 50

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('reference_in', self.option('reference_in').prop['path']))
        self.logger.debug('{} - {}'.format('threads', self.option('threads')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('threads')
        self._memory = '{}G'.format(os.path.getsize(self.option('reference_in').prop['path'])/1024**3+30)

    def end(self):
        super(Bowtie2BuildAgent, self).end()

class Bowtie2BuildTool(Tool):
    def __init__(self, config):
        super(Bowtie2BuildTool, self).__init__(config)
        self.bowtie2_build = 'bioinfo/align/bowtie2-2.2.9/bowtie2-build'

    def run(self):
        super(Bowtie2BuildTool, self).run()
        self.run_bowtie2_build()
        self.set_output()
        self.end()

    def run_bowtie2_build(self):
        cmd = '{} --threads {} {} {}'.format(
            self.bowtie2_build,
            self.option('threads'),
            self.option('reference_in').prop['path'],
            self.option('bt2_index_base')
        )
        cmd_name = 'bowtie2_build'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, ignore_error=True):
        command = self.add_command(cmd_name, cmd, ignore_error=ignore_error)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code in [1, -6, -9]:  # add memory limit by shicaiping at 20190619
            self.add_state("memory_limit", "memory is low!")
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
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bowtie2_build_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.bowtie2_build',
            'instant': False,
            'options': {
                'reference_in': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/Homo_sapiens.GRCh38.94.fa',
                'bt2_index_base': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/Homo_sapiens.GRCh38.94',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()