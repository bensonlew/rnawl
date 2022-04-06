# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class CufflinksAgent(Agent):
    '''
    last_modify: 2019.02.12
    '''
    def __init__(self, parent):
        super(CufflinksAgent, self).__init__(parent)
        options = [
            {'name': 'fr_stranded', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'strand_direct', 'type': 'string', 'default': 'none'},
            {'name': 'sample_name', 'type': 'string', 'default': ''},
            {'name': 'sample_bam', 'type': 'infile', 'format': 'lnc_rna.bam'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'cpu', 'type': 'int', 'default': 8},
            {'name': 'F', 'type': 'float', 'default': 0.1},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'fpkm_cut', 'type': 'int', 'default': 0.1},
            {'name': 'sample_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'sample_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
        ]
        self.add_option(options)
        self.step.add_steps('cufflinks')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.cufflinks.start()
        self.step.update()

    def stepfinish(self):
        self.step.cufflinks.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('fr_stranded', self.option('fr_stranded')))
        if self.option('fr_stranded') == '':
            raise OptionError('whether employ stranded data must be determined')
        self.logger.debug('{} - {}'.format('strand_direct', self.option('strand_direct')))
        if self.option('strand_direct') == '':
            raise OptionError('stranded library must be specified')
        self.logger.debug('{} - {}'.format('sample_name', self.option('sample_name')))
        if self.option('sample_name') == '':
            raise OptionError('sample name must be provided')
        self.logger.debug('{} - {}'.format('sample_bam', self.option('sample_bam').prop['path']))
        if not self.option('sample_bam').is_set:
            raise OptionError('input BAM must be provided')
        self.logger.debug('{} - {}'.format('ref_gtf', self.option('ref_gtf').prop['path']))
        if not self.option('ref_gtf').is_set:
            raise OptionError('reference annotation GFF must be provided')
        self.logger.debug('{} - {}'.format('cpu', self.option('cpu')))
        if self.option('cpu') == '':
            raise OptionError('number of threads must be specified')
        self.logger.debug('{} - {}'.format('F', self.option('F')))
        if self.option('F') == '':
            raise OptionError('minimum isoform fraction must be specified')
        self.logger.debug('{} - {}'.format('fpkm_cut', self.option('fpkm_cut')))
        if self.option('fpkm_cut') == '':
            raise OptionError('FPKM cut-off must be specified')
        self.logger.debug('{} - {}'.format('ref_fa', self.option('ref_fa').prop['path']))
        if not self.option('ref_fa').is_set:
            raise OptionError('genomic seqs FASTA must be provided')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('cpu')
        self._memory = '16G'

    def end(self):
        super(CufflinksAgent, self).end()


class CufflinksTool(Tool):
    def __init__(self, config):
        super(CufflinksTool, self).__init__(config)
        self.cufflinks = 'bioinfo/rna/cufflinks-2.2.1/cufflinks'
        self.sh = 'program/sh'
        self.python = 'miniconda2/bin/python'
        self.gffread = 'bioinfo/rna/cufflinks-2.2.1/gffread'
        self.fasta_range_sh = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/fasta_range.sh')
        self.bioawk_path = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/seq/bioawk/')
        self.filter_gtf_by_range_sh = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_gtf_by_range.sh')
        self.bedtools_path = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/seq/bedtools-2.25.0/bin/')
        self.filter_gtf_by_bed_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_gtf_by_bed.py')
        self.filter_gtf_by_fpkm_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_gtf_by_fpkm.py')
        self.sample_dir = os.path.join(self.work_dir, self.option('sample_name'))
        self.filter_bed = os.path.join(self.work_dir, '{}.filter.bed'.format(self.option('sample_name')))
        self.transcrpts_gtf = os.path.join(self.sample_dir, 'transcripts.gtf')
        self.transcrpts_filter_range_gtf = os.path.join(self.sample_dir, 'transcripts.filter_range.gtf')
        self.transcrpts_filter_fpkm_gtf = os.path.join(self.sample_dir, 'transcripts.filter_fpkm.gtf')
        self.output_fa = os.path.join(self.work_dir, '{}_out.fa'.format(self.option('sample_name')))

    def run(self):
        super(CufflinksTool, self).run()
        self.run_cufflinks()
        self.run_fasta_range()
        self.run_filter_gtf_by_range()
        self.run_filter_gtf_by_bed()
        self.run_filter_gtf_by_fpkm()
        self.run_gffread()
        self.set_output()
        self.end()

    def run_cufflinks(self):
        if self.option('fr_stranded') == 'fr-unstranded':
            cmd = '{} {}'.format(self.cufflinks, self.option('sample_bam').prop['path'])
            cmd += ' -o {}'.format(self.sample_dir)
            cmd += ' -p {}'.format(self.option('cpu'))
            cmd += ' -g {}'.format(self.option('ref_gtf').prop['path'])
            cmd += ' --library-type {}'.format(self.option('fr_stranded'))
            cmd += ' -F {}'.format(self.option('F'))
        elif self.option('fr_stranded') == 'fr-stranded':

            if self.option('strand_direct') == '':
                self.option('strand_direct', 'fr-firststrand')
            else:
                self.option('strand_direct', 'fr-secondstrand')

            cmd = '{} {}'.format(self.cufflinks, self.option('sample_bam').prop['path'])
            cmd += ' -o {}'.format(self.sample_dir)
            cmd += ' -p {}'.format(self.option('cpu'))
            cmd += ' -g {}'.format(self.option('ref_gtf').prop['path'])
            cmd += ' --library-type {}'.format(self.option('strand_direct'))
            cmd += ' -F {}'.format(self.option('F'))
        cmd_name = 'run_cufflinks'
        self.run_code(cmd_name, cmd)

    def run_fasta_range(self):
        cmd = '{} {} {} {} {}'.format(
            self.sh,
            self.fasta_range_sh,
            self.bioawk_path,
            self.option('ref_fa').prop['path'],
            self.filter_bed
        )
        cmd_name = 'run_fasta_range'
        self.run_code(cmd_name, cmd)

    def run_filter_gtf_by_range(self):
        cmd = '{} {} {} {} {} {}'.format(
            self.sh,
            self.filter_gtf_by_range_sh,
            self.bedtools_path,
            self.filter_bed,
            self.transcrpts_gtf,
            self.transcrpts_filter_range_gtf
        )
        cmd_name = 'run_filter_gtf_by_range'
        self.run_code(cmd_name, cmd)

    def run_filter_gtf_by_bed(self):
        cmd = '{} {} {} {} {}'.format(
            self.python,
            self.filter_gtf_by_bed_py,
            self.transcrpts_gtf,
            self.transcrpts_filter_range_gtf,
            self.filter_bed,
        )
        cmd_name = 'run_filter_gtf_by_bed'
        if os.path.getsize(self.transcrpts_filter_range_gtf) == 0:
            self.run_code(cmd_name, cmd)

    def run_filter_gtf_by_fpkm(self):
        cmd = '{} {}'.format(self.python, self.filter_gtf_by_fpkm_py)
        cmd += ' -i {}'.format(self.transcrpts_filter_range_gtf)
        cmd += ' -c {}'.format(self.option('fpkm_cut'))
        cmd += ' -o {}'.format(self.transcrpts_filter_fpkm_gtf)
        cmd_name = 'run_filter_gtf_by_fpkm'
        self.run_code(cmd_name, cmd)

    def run_gffread(self):
        cmd = '{} {} '.format(self.gffread, self.transcrpts_filter_fpkm_gtf)
        cmd += '-g {} '.format(self.option('ref_fa').prop['path'])
        cmd += '-w {}'.format(self.output_fa)
        cmd_name = 'run_gffread'
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
        sample_gtf = os.path.join(self.output_dir, '{}_out.gtf'.format(self.option('sample_name')))
        if os.path.exists(sample_gtf):
            os.remove(sample_gtf)
        os.link(self.transcrpts_filter_fpkm_gtf, sample_gtf)
        self.logger.info('succeed in linking {} to {}'.format(self.transcrpts_filter_fpkm_gtf, sample_gtf))
        self.option('sample_gtf', sample_gtf)
        sample_fa = os.path.join(self.output_dir, '{}_out.fa'.format(self.option('sample_name')))
        if os.path.exists(sample_fa):
            os.remove(sample_fa)
        os.link(self.output_fa, sample_fa)
        self.logger.info('succeed in linking {} to {}'.format(self.output_fa, sample_fa))
        self.option('sample_fa', sample_fa)
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
            'id': 'cufflinks_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.assemble.cufflinks',
            'instant': False,
            'options': {
                'sample_name': 'Xue',
                'sample_bam': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/clean_data/Xue_1115_S4_L002.cufflinks.sort.bam',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
