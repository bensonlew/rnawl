# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import subprocess
import unittest

class MirnaEditBakAgent(Agent):
    '''
    Detect miRNA editing sites
    '''
    def __init__(self, parent):
        super(MirnaEditBakAgent, self).__init__(parent)
        options = [
            {'name': 'species', 'type': 'string', 'default': ''},
            {'name': 'sample', 'type': 'string', 'default': ''},
            {'name': 'filtered_fq', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'index', 'type': 'string', 'default': ''},
            {'name': 'hairpin_fa', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'mature_fa', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'binomial_xls', 'type': 'outfile', 'format': 'small_rna.common'},
        ]
        self.add_option(options)

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))

        self.logger.debug('{} - {}'.format('species', self.option('species')))
        if self.option('species') == '':
            raise OptionError('specie abbreviation must be provided')

        self.logger.debug('{} - {}'.format('sample', self.option('sample')))
        if self.option('sample') == '':
            raise OptionError('sample name must be provided')

        self.logger.debug('{} - {}'.format('filtered_fq', self.option('filtered_fq').prop['path']))
        if not self.option('filtered_fq').is_set:
            raise OptionError('trimmed fastq must be provided')

        self.logger.debug('{} - {}'.format('index', self.option('index')))
        if self.option('index') == '':
            raise OptionError('bowtie index location must be provided')

        self.logger.debug('{} - {}'.format('hairpin_fa', self.option('hairpin_fa').prop['path']))
        if not self.option('hairpin_fa').is_set:
            raise OptionError('pre-miRNA fasta from mirbase must be provided')

        self.logger.debug('{} - {}'.format('mature_fa', self.option('mature_fa').prop['path']))
        if not self.option('mature_fa').is_set:
            raise OptionError('mature-miRNA fasta from mirbase must be provided')

        self.logger.debug('{} - {}'.format('binomial_xls', self.option('binomial_xls')))

        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 4
        self._memory = '24G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', 'mirna_edit_tool_output_dir']
        ])
        super(MirnaEditBakAgent, self).end()

class MirnaEditBakTool(Tool):
    '''
    Detect miRNA editing sites
    '''
    def __init__(self, config):
        super(MirnaEditBakTool, self).__init__(config)
        # set environment variables
        self.anaconda2 = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/miRNA/anaconda2/bin')
        self.set_environ(PATH=self.anaconda2)
        # set program path, '/mnt/ilustre/users/sanger-dev/app' will be attached to the frontier of path
        self.bowtie = 'bioinfo/miRNA/anaconda2/bin/bowtie'
        self.python = 'miniconda2/bin/python'
        self.rnafold = 'bioinfo/miRNA/anaconda2/bin/RNAfold'
        self.perl = 'program/perl-5.24.0/bin/perl'
        # set file path
        self.mir_ref_process = os.path.join(self.config.PACKAGE_DIR, 'small_rna/mirna_ref_process.py')
        self.filtered_hit = os.path.join(self.work_dir, '{}.filtered.hit'.format(self.option('species')))
        self.hairpin_fa = os.path.join(self.work_dir, '{}.hairpin.fa'.format(self.option('species')))
        self.hairpin_rnafold = os.path.join(self.work_dir, '{}.hairpin.rnafold'.format(self.option('species')))
        self.hairpin_u2t_fa = os.path.join(self.work_dir, '{}.hairpin.u2t.fa'.format(self.option('species')))
        self.hairpin_u2t_hit = os.path.join(self.work_dir, '{}.hairpin.u2t.hit'.format(self.option('species')))
        self.mature_fa = os.path.join(self.work_dir, '{}.mature.fa'.format(self.option('species')))
        self.mutation_pl = os.path.join(self.config.PACKAGE_DIR, 'small_rna/analyze_mutation_modified.pl')
        self.mutation_txt = os.path.join(self.work_dir, '{}.mutation.txt'.format(self.option('species')))
        self.binomial_pl = os.path.join(self.config.PACKAGE_DIR, 'small_rna/binomial_analysis_modified.pl')
        self.binomial_xls = os.path.join(self.work_dir, '{}.binomial.xls'.format(self.option('species')))

    def run(self):
        '''
        define running logic
        '''
        super(MirnaEditBakTool, self).run()
        self.filtered_fq_bowtie()
        self.hairpin_filter()
        self.hairpin_fa_rnafold()
        self.hairpin_filter_u2t()
        self.hairpin_u2t_fa_bowtie()
        self.mature_filter()
        self.analyze_mutation()
        self.binomial_analysis()
        self.set_output()
        self.end()

    def filtered_fq_bowtie(self):
        cmd = '{} --trim3 2 --seedmms 1 --maqerr 50 --all -m 1 --best --strata --threads 4 '.format(self.bowtie)
        cmd += '{} '.format(self.option('index'))
        cmd += '{} '.format(self.option('filtered_fq').prop['path'])
        cmd += '{}'.format(self.filtered_hit)
        cmd_name = 'filtered_fq_bowtie'
        self.run_code(cmd_name, cmd)

    def hairpin_filter(self):
        cmd = '{} {} '.format(self.python, self.mir_ref_process)
        cmd += '-i {} '.format(self.option('hairpin_fa').prop['path'])
        cmd += '-s {} '.format(self.option('species'))
        cmd += '-o {}'.format(self.hairpin_fa)
        cmd_name = 'hairpin_filter'
        self.run_code(cmd_name, cmd)

    def hairpin_fa_rnafold(self):
        cmd = '{} --noPS '.format(self.rnafold)
        cmd += '-i {} '.format(self.hairpin_fa)
        cmd += '--outfile={}'.format(os.path.basename(self.hairpin_rnafold))
        cmd_name = 'hairpin_fa_rnafold'
        self.run_code(cmd_name, cmd)

    def hairpin_filter_u2t(self):
        cmd = '{} {} -t '.format(self.python, self.mir_ref_process)
        cmd += '-i {} '.format(self.option('hairpin_fa').prop['path'])
        cmd += '-s {} '.format(self.option('species'))
        cmd += '-o {}'.format(self.hairpin_u2t_fa)
        cmd_name = 'hairpin_filter_u2t'
        self.run_code(cmd_name, cmd)

    def hairpin_u2t_fa_bowtie(self):
        cmd = '{} -f --all -m 1 --best --strata --threads 2 '.format(self.bowtie)
        cmd += '{} '.format(self.option('index'))
        cmd += '{} '.format(self.hairpin_u2t_fa)
        cmd += '{}'.format(self.hairpin_u2t_hit)
        cmd_name = 'hairpin_u2t_fa_bowtie'
        self.run_code(cmd_name, cmd)

    def mature_filter(self):
        cmd = '{} {} '.format(self.python, self.mir_ref_process)
        cmd += '-i {} '.format(self.option('mature_fa').prop['path'])
        cmd += '-s {} '.format(self.option('species'))
        cmd += '-o {}'.format(self.mature_fa)
        cmd_name = 'mature_filter'
        self.run_code(cmd_name, cmd)

    def analyze_mutation(self):
        cmd = '{} {} '.format(self.perl, self.mutation_pl)
        cmd += '{} '.format(self.option('species')) # spc
        cmd += '{} '.format(self.filtered_hit) # spc.filtered.hit
        cmd += '{} '.format(self.hairpin_rnafold) # spc.hairpin.rnafold
        cmd += '{} '.format(self.hairpin_u2t_hit) # spc.hairpin.u2t.hit
        cmd += '{} '.format(self.mature_fa) # spc.mature.fa
        cmd += '{}'.format(self.mutation_txt) # spc.mutation.txt
        cmd_name = 'analyze_mutation'
        self.run_code(cmd_name, cmd)

    def binomial_analysis(self):
        cmd = '{} {} '.format(self.perl, self.binomial_pl)
        cmd += '{} '.format(self.mutation_txt) # spc.mutation.txt
        cmd += '{} '.format(self.hairpin_u2t_hit) # spc.hairpin.u2t.hit
        cmd += '{}'.format(self.binomial_xls) # spc.binomial.xls
        cmd_name = 'binomial_analysis'
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
        '''
        link result in work_dir of tool to output_dir of tool
        '''
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        source = self.binomial_xls
        link_name = os.path.join(self.output_dir, '{}.{}'.format(self.option('sample'), os.path.basename(source)))
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('binomial_xls', link_name)
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_mmu(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'mirna_edit_bak_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'small_rna.mirna_edit_bak',
            'instant': False,
            'options': {
                'species': 'mmu',
                'filtered_fq': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/mirna_edit/mmu.filtered.fq',
                'sample': 'test',
                'index': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean_index',
                'hairpin_fa': '/mnt/ilustre/users/sanger-dev/app/database/mirbase/edit/hairpin.fa',
                'mature_fa': '/mnt/ilustre/users/sanger-dev/app/database/mirbase/edit/mature.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()