# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
import os
from biocluster.tool import Tool
import unittest

class SamtoolsFaidxAgent(Agent):
    def __init__(self, parent):
        super(SamtoolsFaidxAgent, self).__init__(parent)
        options = [
            {'name': 'file_fa', 'type': 'infile', 'format': 'ref_genome_db.fasta'},
            {'name': 'file_fa_fai', 'type': 'outfile', 'format': 'ref_genome_db.common'},
        ]
        self.add_option(options)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('file_fa', self.option('file_fa').prop['path']))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 2
        self._memory = '{}G'.format(os.path.getsize(self.option('file_fa').prop['path'])/1024**3+10)

    def end(self):
        super(SamtoolsFaidxAgent, self).end()

class SamtoolsFaidxTool(Tool):
    def __init__(self, config):
        super(SamtoolsFaidxTool, self).__init__(config)
        self.samtools = 'bioinfo/align/samtools-1.7/samtools'
        self.file_fa_fai = '{}.fai'.format(self.option('file_fa').prop['path'])

    def run(self):
        super(SamtoolsFaidxTool, self).run()
        self.run_samtools_faidx()
        self.set_output()
        self.end()

    def run_samtools_faidx(self):
        cmd = '{} faidx {}'.format(
            self.samtools,
            self.option('file_fa').prop['path']
        )
        cmd_name = 'samtools_faidx'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
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
        source = self.file_fa_fai
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('file_fa_fai', link_name)
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
            'id': 'samtool_faidx_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_genome_db.samtools_faidx',
            'instant': False,
            'options': {
                'file_fa': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_genome_db/Homo_sapiens.GRCh38.94.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()