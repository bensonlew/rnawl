# -*- coding: utf-8 -*-
# __author__ = 'litangjian, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class SnptoolsAgent(Agent):
    '''
    last_modify: 2019.02.28
    '''
    def __init__(self, parent):
        super(SnptoolsAgent, self).__init__(parent)
        options = [
            {'name': 'bam_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'adjust_MQ', 'type': 'int', 'default': 50},
            {'name': 'ref_fa', 'type': 'infile', 'format':'lnc_rna.fasta'},
            {'name': 'chr_bed', 'type': 'infile', 'format':'lnc_rna.common'},
            {'name': 'min_MQ', 'type': 'int', 'default': 10},
            {'name': 'min_BQ', 'type': 'int', 'default': 30},
            {'name': 'output_tags', 'type': 'string', 'default': 'DP,AD'},
            {'name': 'output_type', 'type': 'string', 'default': 'v'},
            {'name': 'format_fields', 'type': 'string', 'default': 'GQ,GP'},
            {'name': 'min_QUAL', 'type': 'int', 'default': 10},
            {'name': 'min_avg_DP', 'type': 'int', 'default': 5},
            {'name': 'snp_gap', 'type': 'int', 'default': 3},
            {'name': 'chr_vcf', 'type': 'outfile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('snptools')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.snptools.start()
        self.step.update()

    def step_end(self):
        self.step.snptools.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('bam_list').is_set:
            self.logger.debug('{} = {}'.format('bam_list', self.option('bam_list').prop['path']))
        else:
            raise OptionError('BAM list file must be provided')
        self.logger.debug('{} = {}'.format('adjust_MQ', self.option('adjust_MQ')))
        if self.option('ref_fa').is_set:
            self.logger.debug('{} = {}'.format('ref_fa', self.option('ref_fa').prop['path']))
            self.infile_size = os.path.getsize(self.option('ref_fa').prop['path'])
        else:
            raise OptionError('reference FASTA file must be provided')
        if self.option('chr_bed').is_set:
            self.logger.debug('{} = {}'.format('chr_bed', self.option('chr_bed').prop['path']))
        else:
            raise OptionError('positions included BED file must be provided')
        self.logger.debug('{} = {}'.format('min_MQ', self.option('min_MQ')))
        self.logger.debug('{} = {}'.format('min_BQ', self.option('min_BQ')))
        self.logger.debug('{} = {}'.format('output_tags', self.option('output_tags')))
        self.logger.debug('{} = {}'.format('output_type', self.option('output_type')))
        self.logger.debug('{} = {}'.format('format_fields', self.option('format_fields')))
        self.logger.debug('{} = {}'.format('min_QUAL', self.option('min_QUAL')))
        self.logger.debug('{} = {}'.format('min_avg_DP', self.option('min_avg_DP')))
        self.logger.debug('{} = {}'.format('snp_gap', self.option('snp_gap')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 4 + 4))

    def end(self):
        super(SnptoolsAgent, self).end()

class SnptoolsTool(Tool):
    def __init__(self, config):
        super(SnptoolsTool, self).__init__(config)
        self.samtools = 'bioinfo/align/samtools-1.6/samtools'
        self.bcftools = 'bioinfo/align/bcftools-1.6/bcftools'
        self.chr = os.path.basename(self.option('chr_bed').prop['path']).replace('.bed', '')
        self.chr_bcf = os.path.join(self.work_dir, '{}.bcf'.format(self.chr))
        self.chr_call_vcf = os.path.join(self.work_dir, '{}.call.vcf'.format(self.chr))
        self.chr_filter_vcf = os.path.join(self.work_dir, '{}.filter.vcf'.format(self.chr))

    def run(self):
        super(SnptoolsTool, self).run()
        self.run_samtools_mpileup()
        self.run_bcftools_call()
        self.run_bcftools_filter()
        self.set_output()
        self.end()

    def run_samtools_mpileup(self):
        cmd = '{} mpileup -g -u'.format(self.samtools)
        cmd += ' -b {}'.format(self.option('bam_list').prop['path'])
        cmd += ' -C {}'.format(self.option('adjust_MQ'))
        cmd += ' -f {}'.format(self.option('ref_fa').prop['path'])
        cmd += ' -l {}'.format(self.option('chr_bed').prop['path'])
        cmd += ' -q {}'.format(self.option('min_MQ'))
        cmd += ' -Q {}'.format(self.option('min_BQ'))
        cmd += ' -o {}'.format(self.chr_bcf)
        cmd += ' -t {}'.format(self.option('output_tags'))
        cmd_name = 'run_samtools_mpileup'
        self.run_code(cmd_name, cmd)

    def run_bcftools_call(self):
        cmd = '{} call -v -m {}'.format(self.bcftools, self.chr_bcf)
        cmd += ' -o {}'.format(self.chr_call_vcf)
        cmd += ' -O {}'.format(self.option('output_type'))
        cmd += ' -f {}'.format(self.option('format_fields'))
        cmd_name = 'run_bcftools_call'
        self.run_code(cmd_name, cmd)

    def run_bcftools_filter(self):
        cmd = '{} filter {}'.format(self.bcftools, self.chr_call_vcf)
        cmd += " -e 'QUAL<{} || DP<{}'".format(
            self.option('min_QUAL'), self.option('min_avg_DP') * len(list(open(self.option('bam_list').prop['path'])))
        )
        cmd += ' -g {}'.format(self.option('snp_gap'))
        cmd += ' -o {}'.format(self.chr_filter_vcf)
        cmd += ' -O {}'.format(self.option('output_type'))
        cmd_name = 'run_bcftools_filter'
        self.run_code(cmd_name, cmd, shell=True)

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
        chr_vcf = os.path.join(self.output_dir, os.path.basename(self.chr_filter_vcf))
        if os.path.exists(chr_vcf):
            os.remove(chr_vcf)
        os.link(self.chr_filter_vcf, chr_vcf)
        self.logger.info('succeed in linking {} to {}'.format(self.chr_filter_vcf, chr_vcf))
        self.option('chr_vcf').set_path(chr_vcf)
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
            'id': 'snptools_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.snptools',
            'instant': False,
            'options': {
                'bam_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/test_data/bam_list.txt',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
                'chr_bed': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/split_bed/chr1.bed',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()