# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import unittest

class SalmonAgent(Agent):
    '''
    last_modify: 2019.04.09
    '''
    def __init__(self, parent):
        super(SalmonAgent, self).__init__(parent)
        options = [
            {'name': 'fq_type', 'type': 'string', 'default': None},  # PE SE
            {'name': 'transcripts', 'type': 'infile', 'format': 'ref_rna_v2.fasta'},
            {'name': 'lib_type', 'type': 'string', 'default': 'A'},
            {'name': 'mates1', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'mates2', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'unmated_reads', 'type': 'infile', 'format': 'ref_rna_v2.fastq'},
            {'name': 'threads', 'type': 'int', 'default': 10},
            {'name': 't2g', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'quant_sf', 'type': 'outfile', 'format': 'lnc_rna.common'},
            {'name': 'quant_genes_sf', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('salmon')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.salmon.start()
        self.step.update()

    def stepfinish(self):
        self.step.salmon.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('threads')
        if self.option('fq_type') == 'PE':
            size = os.path.getsize(self.option('mates1').path) + os.path.getsize(self.option('mates2').path)
        elif self.option('fq_type') == 'SE':
            size = os.path.getsize(self.option('unmated_reads').path)
        self._memory = '{}G'.format(int(float(size) / 1024 ** 3 + 10))

    def end(self):
        super(SalmonAgent, self).end()

class SalmonTool(Tool):
    def __init__(self, config):
        super(SalmonTool, self).__init__(config)
        self.salmon = 'bioinfo/rna/Salmon-0.8.2_linux_x86_64/bin/salmon'
        self.index = os.path.join(self.work_dir, 'index')
        self.output = os.path.join(self.work_dir, '{}'.format(self.option('sample_name')))
        self.quant_sf = os.path.join(self.output, 'quant.sf')
        self.quant_genes_sf = os.path.join(self.output, 'quant.genes.sf')

    def run(self):
        super(SalmonTool, self).run()
        self.run_salmon_index()
        if self.option('fq_type') == 'PE':
            self.run_salmon_quant_pe()
        elif self.option('fq_type') == 'SE':
            self.run_salmon_quant_se()
        self.set_output()
        self.end()

    def run_salmon_index(self):
        cmd = '{} index'.format(self.salmon)
        cmd += ' -t {}'.format(self.option('transcripts').path)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -p {}'.format(self.option('threads'))
        cmd_name = 'run_salmon_index'
        self.run_code(cmd_name, cmd)

    def run_salmon_quant_pe(self):
        cmd = '{} quant --gcBias'.format(self.salmon)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -l {}'.format(self.option('lib_type'))
        cmd += ' -1 {}'.format(self.option('mates1').path)
        cmd += ' -2 {}'.format(self.option('mates2').path)
        cmd += ' -o {}'.format(self.output)
        cmd += ' -p {}'.format(self.option('threads'))
        cmd += ' -g {}'.format(self.option('t2g').path)
        cmd_name = 'run_salmon_quant_pe'
        self.run_code(cmd_name, cmd)

    def run_salmon_quant_se(self):
        cmd = '{} quant --gcBias'.format(self.salmon)
        cmd += ' -i {}'.format(self.index)
        cmd += ' -l {}'.format(self.option('lib_type'))
        cmd += ' -r {}'.format(self.option('unmated_reads').path)
        cmd += ' -o {}'.format(self.output)
        cmd += ' -p {}'.format(self.option('threads'))
        cmd += ' -g {}'.format(self.option('t2g').path)
        cmd_name = 'run_salmon_quant_se'
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
        source = self.quant_sf
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('quant_sf').set_path(link_name)
        source = self.quant_genes_sf
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('quant_genes_sf').set_path(link_name)
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
            'id': 'salmon_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.salmon',
            'instant': False,
            'options': {
                'fq_type': 'PE',
                'transcripts': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'mates1': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_l.fastq',
                'mates2': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_r.fastq',
                't2g': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/t2g.pairs',
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
            'id': 'salmon_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.expression.salmon',
            'instant': False,
            'options': {
                'fq_type': 'SE',
                'transcripts': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/all_transcripts.fa',
                'unmated_reads': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/sickle_dir/Ctr_Liver_1_sickle_l.fastq',
                't2g': '/mnt/lustre/users/sanger/sg-users/qinjincheng/ref_rna_v2/expression/NewTranscripts/t2g.pairs',
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
