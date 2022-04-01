# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class KallistoAgent(Agent):
    '''
    last_modify: 2019.04.09
    '''
    def __init__(self, parent):
        super(KallistoAgent, self).__init__(parent)
        options = [
            {'name': 'fq_type', 'type': 'string', 'default': None}, # PE SE
            {'name': 'fasta_file', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'threads', 'type': 'int', 'default': 10},
            {'name': 'fastq_l_file', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'fastq_r_file', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'fastq_file', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'strand_specific', 'type': 'bool', 'default': False},
            {'name': 'stranded', 'type': 'string', 'default': None}, # rf (workflow forward) fr (workflow reverse)
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'abundance_tsv', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'kmer', 'type': 'int', 'default': 31}
        ]
        self.add_option(options)
        self.step.add_steps('kallisto')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kallisto.start()
        self.step.update()

    def stepfinish(self):
        self.step.kallisto.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('threads')
        if self.option('fq_type') == 'PE':
            size = os.path.getsize(self.option('fastq_l_file').path) + os.path.getsize(self.option('fastq_r_file').path)
        elif self.option('fq_type') == 'SE':
            size = os.path.getsize(self.option('fastq_file').path)
        self._memory = '{}G'.format(int(float(size) / 1024 ** 3 + 10))

    def end(self):
        super(KallistoAgent, self).end()

class KallistoTool(Tool):
    def __init__(self, config):
        super(KallistoTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.kallisto = 'bioinfo/rna/kallisto_linux-v0.43.1/kallisto'
        self.index = os.path.join(self.work_dir, 'index')
        self.output = os.path.join(self.work_dir, self.option('sample_name'))
        self.abundance_tsv = os.path.join(self.output, 'abundance.tsv')

    def run(self):
        super(KallistoTool, self).run()
        self.run_kallisto_index()
        if self.option('fq_type') == 'PE':
            self.run_kallisto_quant_pe()
        elif self.option('fq_type') == 'SE':
            self.run_kallisto_quant_se()
        self.set_output()
        self.end()

    def run_kallisto_index(self):
        cmd = '{} index'.format(self.kallisto)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -k {}'.format(self.option('kmer'))
        cmd += ' {}'.format(self.option('fasta_file').path)
        cmd_name = 'run_kallisto_index'
        self.run_code(cmd_name, cmd)

    def run_kallisto_quant_pe(self):
        cmd = '{} quant'.format(self.kallisto)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -o {}'.format(self.output)
        if self.option('strand_specific'):
            cmd += ' --{}-stranded'.format(self.option('stranded'))
        cmd += ' -t {}'.format(self.option('threads'))
        cmd += ' {} {}'.format(self.option('fastq_l_file').path, self.option('fastq_r_file').path)
        cmd_name = 'run_kallisto_quant_pe'
        self.run_code(cmd_name, cmd)

    def run_kallisto_quant_se(self):
        if self.option('fastq_file').prepare_kallisto():
            cmd = '{} quant --single'.format(self.kallisto)
            cmd += ' -i {}'.format(self.index)
            cmd += ' -o {}'.format(self.output)
            cmd += ' -l {}'.format(self.option('fastq_file').kallisto['l'])
            cmd += ' -s {}'.format(self.option('fastq_file').kallisto['s'])
            cmd += ' -t {}'.format(self.option('threads'))
            cmd += ' {}'.format(self.option('fastq_file').path)
            cmd_name = 'run_kallisto_quant_se'
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
        source = self.abundance_tsv
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('abundance_tsv').set_path(link_name)
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
            'id': 'kallisto_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.kallisto',
            'instant': False,
            'options': {
                'fq_type': 'PE',
                'fasta_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'fastq_l_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_l.fastq',
                'fastq_r_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_r.fastq',
                'sample_name': 'Ctr_Liver_1'
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
            'id': 'kallisto_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.kallisto',
            'instant': False,
            'options': {
                'fq_type': 'SE',
                'fasta_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'fastq_file': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_l.fastq',
                'sample_name': 'Ctr_Liver_1'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_se')])
    unittest.TextTestRunner(verbosity=2).run(suite)
