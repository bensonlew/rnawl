# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class TrimmomaticAgent(Agent):
    '''
    last_modify: 2019.04.29
    '''
    def __init__(self, parent):
        super(TrimmomaticAgent, self).__init__(parent)
        options = [
            {'name': 'fq_type', 'type': 'string', 'default': None}, # PE SE
            {'name': 'threads', 'type': 'int', 'default': 8},
            {'name': 'phred', 'type': 'int', 'default': 33},
            {'name': 'input_file1', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'input_file2', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'input_file', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'seed_mismatches', 'type': 'int', 'default': 2},
            {'name': 'palindrome_clip_threshold', 'type': 'int', 'default': 30},
            {'name': 'simple_clip_threshold', 'type': 'int', 'default': 10},
            {'name': 'leading_quality', 'type': 'int', 'default': 3},
            {'name': 'trailing_quality', 'type': 'int', 'default': 3},
            {'name': 'window_size', 'type': 'int', 'default': 4},
            {'name': 'required_quality', 'type': 'int', 'default': 15},
            {'name': 'minimum_length', 'type': 'int', 'default': 36},
            {'name': 'output_file_1p', 'type': 'outfile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'output_file_2p', 'type': 'outfile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'output_file', 'type': 'outfile', 'format': 'ref_rna_v2.fastq'},
        ]
        self.add_option(options)
        self.step.add_steps('trimmomatic')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.trimmomatic.start()
        self.step.update()

    def stepfinish(self):
        self.step.trimmomatic.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('threads')
        if self.option('fq_type') == 'PE':
            size = os.path.getsize(self.option('input_file1').path) + os.path.getsize(self.option('input_file2').path)
        elif self.option('fq_type') == 'SE':
            size = os.path.getsize(self.option('input_file').path)
        self._memory = '{}G'.format(int(float(size) / 1024 ** 4 + 8))

    def end(self):
        super(TrimmomaticAgent, self).end()

class TrimmomaticTool(Tool):
    def __init__(self, config):
        super(TrimmomaticTool, self).__init__(config)
        self.trimmomatic = 'bioinfo/rna/miniconda2/share/trimmomatic-0.39-1/trimmomatic'
        if self.option('fq_type') == 'PE':
            self.output_file_1p = os.path.join(self.output_dir, '{}_1P.fq.gz'.format(self.option('sample_name')))
            self.output_file_1u = os.path.join(self.output_dir, '{}_1U.fq.gz'.format(self.option('sample_name')))
            self.output_file_2p = os.path.join(self.output_dir, '{}_2P.fq.gz'.format(self.option('sample_name')))
            self.output_file_2u = os.path.join(self.output_dir, '{}_2U.fq.gz'.format(self.option('sample_name')))
        elif self.option('fq_type') == 'SE':
            self.output_file = os.path.join(self.output_dir, '{}.fq.gz'.format(self.option('sample_name')))
        self.adapters_fa = os.path.join(
            self.config.SOFTWARE_DIR,
            'bioinfo/rna/miniconda2/share/trimmomatic-0.39-1/adapters/Adapters.fa'
        )

    def run(self):
        super(TrimmomaticTool, self).run()
        if self.option('fq_type') == 'PE':
            self.run_trimmomatic_pe()
        elif self.option('fq_type') == 'SE':
            self.run_trimmomatic_se()
        self.set_output()
        self.end()

    def run_trimmomatic_pe(self):
        cmd = '{} PE'.format(self.trimmomatic)
        cmd += ' -threads {}'.format(self.option('threads'))
        cmd += ' -phred{}'.format(self.option('phred'))
        cmd += ' {}'.format(self.option('input_file1').path)
        cmd += ' {}'.format(self.option('input_file2').path)
        cmd += ' {}'.format(self.output_file_1p)
        cmd += ' {}'.format(self.output_file_1u)
        cmd += ' {}'.format(self.output_file_2p)
        cmd += ' {}'.format(self.output_file_2u)
        cmd += ' ILLUMINACLIP:{}:{}:{}:{}'.format(
            self.adapters_fa,
            self.option('seed_mismatches'),
            self.option('palindrome_clip_threshold'),
            self.option('simple_clip_threshold')
        )
        cmd += ' LEADING:{}'.format(self.option('leading_quality'))
        cmd += ' TRAILING:{}'.format(self.option('trailing_quality'))
        cmd += ' SLIDINGWINDOW:{}:{}'.format(self.option('window_size'), self.option('required_quality'))
        cmd += ' MINLEN:{}'.format(self.option('minimum_length'))
        cmd_name = 'run_trimmomatic_pe'
        self.run_code(cmd_name, cmd)

    def run_trimmomatic_se(self):
        cmd = '{} SE'.format(self.trimmomatic)
        cmd += ' -threads {}'.format(self.option('threads'))
        cmd += ' -phred{}'.format(self.option('phred'))
        cmd += ' {}'.format(self.option('input_file').path)
        cmd += ' {}'.format(self.output_file)
        cmd += ' ILLUMINACLIP:{}:{}:{}:{}'.format(
            self.adapters_fa,
            self.option('seed_mismatches'),
            self.option('palindrome_clip_threshold'),
            self.option('simple_clip_threshold')
        )
        cmd += ' LEADING:{}'.format(self.option('leading_quality'))
        cmd += ' TRAILING:{}'.format(self.option('trailing_quality'))
        cmd += ' SLIDINGWINDOW:{}:{}'.format(self.option('window_size'), self.option('required_quality'))
        cmd += ' MINLEN:{}'.format(self.option('minimum_length'))
        cmd_name = 'run_trimmomatic_se'
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
                        self.set_error('fail to run %s, abord', variables=(n), code="33712002")

    def set_output(self):
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        if self.option('fq_type') == 'PE':
            self.option('output_file_1p').set_path(self.output_file_1p)
            self.option('output_file_2p').set_path(self.output_file_2p)
        elif self.option('fq_type') == 'SE':
            self.option('output_file').set_path(self.output_file)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_pe(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'trimmomatic_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.trimmomatic',
            'instant': False,
            'options': {
                'fq_type': 'PE',
                'input_file1': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_demo/remote_input/fastq_dir/TR1-3_S6_L007_R1_001.fastq.gz',
                'input_file2': '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_demo/remote_input/fastq_dir/TR1-3_S6_L007_R2_001.fastq.gz',
                'sample_name': 'TR1_3'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_se(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'trimmomatic_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.trimmomatic',
            'instant': False,
            'options': {
                'fq_type': 'SE',
                'input_file': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/chip/00.raw_data/MESC_H3K4me1.fq.gz',
                'sample_name': 'MESC_H3K4me1',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_pe')])
    unittest.TextTestRunner(verbosity=2).run(suite)
