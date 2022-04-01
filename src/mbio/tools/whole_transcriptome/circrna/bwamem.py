# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class BwamemAgent(Agent):
    '''
    last_modify: 2019.8.22
    '''

    def __init__(self, parent):
        super(BwamemAgent, self).__init__(parent)
        options = [
            {'name': 'sample', 'type': 'string', 'default': None},
            {'name': 'genome', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'fq1', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'fq2', 'type': 'infile', 'format': 'whole_transcriptome.fastq'},
            {'name': 'threads', 'type': 'int', 'default': 20},
            {'name': 'rnasam', 'type': 'outfile', 'format': 'whole_transcriptome.common'},

        ]
        self.add_option(options)
        self._memory_increase_step = 100

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 21
        self._memory = '15G'

    def end(self):
        super(BwamemAgent, self).end()


class BwamemTool(Tool):
    def __init__(self, config):
        super(BwamemTool, self).__init__(config)

        self.program = {
            'bwa': 'install_packages/bwa-0.7.16a/bwa',

        }
        self.file = {
            'rnasam': os.path.join(self.output_dir, '{}.sam'.format(self.option('sample')))

        }

    def run(self):
        super(BwamemTool, self).run()
        if self.check_no_index():
            self.run_bwa_index()
        self.run_bwa_mem()
        self.set_output()
        self.end()

    def check_no_index(self):
        for suffix in ('amb', 'ann', 'bwt', 'pac', 'sa'):
            if not os.path.isfile('{}.{}'.format(self.option('genome').path, suffix)):
                return True
        else:
            return False

    def run_bwa_index(self):
        cmd = '{} {} {} {} {}'.format(self.program['bwa'], 'index', '-a', 'bwtsw', self.option('genome').path)
        cmd_name = 'run_bwa_index'
        runcmd(self, cmd_name, cmd)

    def run_bwa_mem(self):
        if self.option('fq2').is_set:
            cmd = '{} {} {} {} {} {} {} {} {} {} {}'.format(self.program['bwa'], 'mem', '-t', self.option('threads'), '-T', 19, self.option('genome').path,
                                                            self.option('fq1').path, self.option('fq2').path, '-o',
                                                            self.file['rnasam'])
            # cmd += ' -t {}'.format(self.option('threads'))
            cmd_name = 'run_bwa_mem'
            if os.path.exists(self.file['rnasam']):
                os.remove(self.file['rnasam'])
                # pass
            # else:
            command = self.add_command(cmd_name, cmd, ignore_error=True)
            command.run()
            self.wait()
            if command.return_code == 0:
                # self.samtools_sort()
                pass
            elif command.return_code == 1:
                self.logger.info("return code: {}".format(command.return_code))
                self.add_state('memory_limit', 'memory is low!')
            elif command.return_code is None:
                command.rerun()
                if command.return_code == 0:
                    self.logger.info("succeed in running run_bwa_mem")
                    pass
            else:
                self.set_error("fail to run run_bwa_mem")
            # runcmd(self, cmd_name, cmd, ignore_error=True)
        else:
            cmd = '{} {} {} {} {} {} {} {} {}'.format(self.program['bwa'], 'mem', '-t', self.option('threads'), '-T', 19, self.option('genome').path,
                                                self.option('fq1').path, '-o',
                                                self.file['rnasam'])
            cmd += ' -t {}'.format(self.option('threads'))
            cmd_name = 'run_bwa_mem'
            if os.path.exists(self.file['rnasam']):
                os.remove(self.file['rnasam'])
                # pass
            command = self.add_command(cmd_name, cmd, ignore_error=True)
            command.run()
            self.wait()
            if command.return_code == 0:
                # self.samtools_sort()
                pass
            elif command.return_code == 1:
                self.logger.info("return code: {}".format(command.return_code))
                self.add_state('memory_limit', 'memory is low!')
            elif command.return_code is None:
                command.rerun()
                if command.return_code == 0:
                    self.logger.info("succeed in running run_bwa_mem")
                    pass
            else:
                self.set_error("fail to run run_bwa_mem")

            # runcmd(self, cmd_name, cmd, ignore_error=True)

    def set_output(self):
        self.option('rnasam').set_path(self.file['rnasam'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bwamem_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.circrna.bwamem',
            'instant': False,
            'options': {
                'genome': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa',
                'fq1': '/mnt/ilustre/users/sanger-dev/workspace/20191010/Longrna_workflow_2617_9288/FastpRna/output/fastq/S1.clean.1.fastq',
                # 'fq2': '/mnt/ilustre/users/sanger-dev/workspace/20191010/Longrna_workflow_2617_9288/FastpRna/output/fastq/S1.clean.2.fastq',
                'sample': 'S1'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
