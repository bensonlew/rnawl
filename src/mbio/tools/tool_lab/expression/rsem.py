# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class RsemAgent(Agent):
    '''
    last_modify: 2019.04.09
    '''
    def __init__(self, parent):
        super(RsemAgent, self).__init__(parent)
        options = [
            {'name': 'fq_type', 'type': 'string', 'default': None}, # PE SE
            {'name': 'upstream_read', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'downstream_read', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'g2t', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'num_threads', 'type': 'int', 'default': 10},
            {'name': 'reference_fasta', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'phred_quals', 'type': 'int', 'default': None}, # 33 64
            # 0 for a strand-specific protocol where all (upstream) reads are derived from the reverse strand
            # 1 for a strand-specific protocol where all (upstream) reads are derived from the forward strand
            # 0.5 for a non-strand-specific protocol
            {'name': 'forward_prob', 'type': 'float', 'default': 0.5}, # 0 (workflow forward) 1 (workflow reverse)
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'isoforms_results', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'genes_results', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'bowtie2-mismatch-rate', 'type': 'float', 'default': 0.1},
            {'name': 'bowtie2-k', 'type': 'int', 'default': 200}
        ]
        self.add_option(options)
        self.step.add_steps('rsem')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.rsem.start()
        self.step.update()

    def stepfinish(self):
        self.step.rsem.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('num_threads')
        size = os.path.getsize(self.option('upstream_read').path)
        if self.option('fq_type') == 'PE':
            size += os.path.getsize(self.option('downstream_read').path)
        self._memory = '{}G'.format(int(float(size) / 1024 ** 3 + 10))

    def end(self):
        super(RsemAgent, self).end()

class RsemTool(Tool):
    def __init__(self, config):
        super(RsemTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        if "dev" in self.config.PACKAGE_DIR:
            self.rsem_prepare_reference = 'bioinfo/rna/RSEM-1.2.31/bin/rsem-prepare-reference'
        else:
            self.rsem_prepare_reference = '/bioinfo/rna/RSEM-1.3.1/rsem-prepare-reference'
        # self.rsem_prepare_reference = 'bioinfo/rna/RSEM-1.2.31/bin/rsem-prepare-reference'
        # self.rsem_prepare_reference = 'bioinfo/rna/RSEM-1.3.1/bin/rsem-prepare-reference'
        self.bowtie2_path = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/align/bowtie2-2.2.9')
        # self.rsem_calculate_expression = 'bioinfo/rna/RSEM-1.2.31/bin/rsem-calculate-expression'
        # self.rsem_calculate_expression = 'bioinfo/rna/RSEM-1.3.1/bin/rsem-calculate-expression'
        if "dev" in self.config.PACKAGE_DIR:
            self.rsem_calculate_expression = '/bioinfo/rna/RSEM-1.2.31/bin/rsem-calculate-expression'
        else:
            self.rsem_calculate_expression = '/bioinfo/rna/RSEM-1.3.1/rsem-calculate-expression'
        self.perl = software_dir + '/program/perl/perls/perl-5.24.0/bin/'
        self.set_environ(PATH=self.perl)
        if "dev" in self.config.PACKAGE_DIR:
            self.rsem_path = software_dir + '/bioinfo/rna/RSEM-1.2.31/bin'
        else:
            self.rsem_path = software_dir + '/bioinfo/rna/RSEM-1.3.1/'
        self.set_environ(PATH=self.rsem_path)
        self.index = os.path.join(self.work_dir, 'index')
        self.sample_name = os.path.join(self.work_dir, self.option('sample_name'))
        self.genes_results = '{}.genes.results'.format(self.sample_name)
        self.isoforms_results = '{}.isoforms.results'.format(self.sample_name)

    def run(self):
        super(RsemTool, self).run()
        self.run_rsem_prepare_reference()
        if self.option('fq_type') == 'PE':
            self.run_rsem_calculate_expression_pe()
        elif self.option('fq_type') == 'SE':
            self.run_rsem_calculate_expression_se()
        self.set_output()
        self.end()

    def run_rsem_prepare_reference(self):
        cmd = '{}'.format(self.rsem_prepare_reference)
        cmd += ' --transcript-to-gene-map {}'.format(self.option('g2t').path)
        cmd += ' --bowtie2 --bowtie2-path {}'.format(self.bowtie2_path)
        cmd += ' --num-threads {}'.format(self.option('num_threads'))
        cmd += ' {} {}'.format(self.option('reference_fasta').path, self.index)
        cmd_name = 'run_rsem_prepare_reference'
        self.run_code(cmd_name, cmd)

    def run_rsem_calculate_expression_pe(self):
        cmd = '{} --paired-end --estimate-rspd'.format(self.rsem_calculate_expression)
        cmd += ' --phred{}-quals'.format(self.option('phred_quals'))
        cmd += ' --bowtie2 --bowtie2-path {}'.format(self.bowtie2_path)
        cmd += ' --num-threads {}'.format(self.option('num_threads'))
        cmd += ' --forward-prob {}'.format(self.option('forward_prob'))
        cmd += ' --bowtie2-mismatch-rate {}'.format(self.option('bowtie2-mismatch-rate'))
        cmd += ' --bowtie2-k {}'.format(self.option('bowtie2-k'))
        cmd += ' {} {}'.format(self.option('upstream_read').path, self.option('downstream_read').path)
        cmd += ' {} {}'.format(self.index, self.sample_name)
        cmd_name = 'run_rsem_calculate_expression_pe'
        self.run_code(cmd_name, cmd)

    def run_rsem_calculate_expression_se(self):
        cmd = '{} --estimate-rspd'.format(self.rsem_calculate_expression)
        cmd += ' --phred{}-quals'.format(self.option('phred_quals'))
        cmd += ' --bowtie2 --bowtie2-path {}'.format(self.bowtie2_path)
        cmd += ' --bowtie2-mismatch-rate {}'.format(self.option('bowtie2-mismatch-rate'))
        cmd += ' --bowtie2-k {}'.format(self.option('bowtie2-k'))
        cmd += ' --num-threads {}'.format(self.option('num_threads'))
        cmd += ' {} {} {}'.format(self.option('upstream_read').path, self.index, self.sample_name)
        cmd_name = 'run_rsem_calculate_expression_se'
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
        source = self.isoforms_results
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('isoforms_results').set_path(link_name)
        source = self.genes_results
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('genes_results').set_path(link_name)
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
            'id': 'rsem_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.rsem',
            'instant': False,
            'options': {
                'fq_type': 'PE',
                'upstream_read': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_l.fastq',
                'downstream_read': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_r.fastq',
                'g2t': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/g2t.pairs',
                'reference_fasta': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'phred_quals': 33,
                'forward_prob': 0.5,
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
            'id': 'rsem_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.rsem',
            'instant': False,
            'options': {
                'fq_type': 'SE',
                'upstream_read': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_l.fastq',
                'g2t': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/g2t.pairs',
                'reference_fasta': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'sample_name': 'Ctr_Liver_1',
                'phred_quals': 33,
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_se')])
    unittest.TextTestRunner(verbosity=2).run(suite)
