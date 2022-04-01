# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import subprocess
import unittest

class MirnaEditAgent(Agent):
    '''
    Detect miRNA editing sites
    '''
    def __init__(self, parent):
        super(MirnaEditAgent, self).__init__(parent)
        options = [
            {'name': 'species', 'type': 'string', 'default': ''},
            {'name': 'sample', 'type': 'string', 'default': ''},
            {'name': 'filtered_fq', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'index', 'type': 'string', 'default': ''},
            {'name': 'hairpin_fa', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'mature_fa', 'type': 'infile', 'format': 'small_rna.common'},
            {'name': 'result_xls', 'type': 'outfile', 'format': 'small_rna.common'},
        ]
        self.add_option(options)
        self._memory_increase_step = 30

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

        self.logger.debug('{} - {}'.format('result_xls', self.option('result_xls')))

        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 4
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', 'mirna_edit_tool_output_dir']
        ])
        super(MirnaEditAgent, self).end()

class MirnaEditTool(Tool):
    '''
    Detect miRNA editing sites
    '''
    def __init__(self, config):
        super(MirnaEditTool, self).__init__(config)
        # set environment variables
        self.anaconda2 = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/miRNA/anaconda2/bin')
        self.set_environ(PATH=self.anaconda2)
        # set program path, '/mnt/ilustre/users/sanger-dev/app' will be attached to the frontier of path
        self.bowtie =  'bioinfo/align/bowtie-1.2.3-linux-x86_64/bowtie'
        # self.bowtie = 'bioinfo/miRNA/anaconda2/bin/bowtie'
        self.python = 'program/Python/bin/python'
        self.perl = 'program/perl-5.24.0/bin/perl'
        # set file path
        self.mirna_ref_process = os.path.join(self.config.PACKAGE_DIR, 'small_rna/mirna_ref_process.py')
        self.mirna_site_filter = os.path.join(self.config.PACKAGE_DIR, 'small_rna/mirna_site_filter.py')
        self.filtered_hit = os.path.join(self.work_dir, '{}.filtered.hit'.format(self.option('species')))
        self.hairpin_fa = os.path.join(self.work_dir, '{}.hairpin.fa'.format(self.option('species')))
        self.hairpin_u2t_fa = os.path.join(self.work_dir, '{}.hairpin.u2t.fa'.format(self.option('species')))
        self.hairpin_u2t_hit = os.path.join(self.work_dir, '{}.hairpin.u2t.hit'.format(self.option('species')))
        self.mature_fa = os.path.join(self.work_dir, '{}.mature.fa'.format(self.option('species')))
        self.analysis_script = os.path.join(self.config.PACKAGE_DIR, 'small_rna/mirna_edit_analysis.py')
        self.mirna_site_json = os.path.join(self.config.SOFTWARE_DIR, 'database/mirbase/miRNA.site.json')
        self.filter_site_json = os.path.join(self.work_dir, 'filter.site.json')
        self.mirna_name_json = os.path.join(self.config.SOFTWARE_DIR, 'database/mirbase/miRNA.name.json')
        self.result_xls = os.path.join(self.work_dir, '{}.result.xls'.format(self.option('species')))

    def run(self):
        '''
        define running logic
        '''
        super(MirnaEditTool, self).run()
        self.filtered_fq_bowtie()
        self.hairpin_filter()
        self.hairpin_filter_u2t()
        self.hairpin_u2t_fa_bowtie()
        self.mature_filter()
        self.mirna_edit_analysis()
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
        cmd = '{} {} '.format(self.python, self.mirna_ref_process)
        cmd += '-i {} '.format(self.option('hairpin_fa').prop['path'])
        cmd += '-s {} '.format(self.option('species'))
        cmd += '-o {}'.format(self.hairpin_fa)
        cmd_name = 'hairpin_filter'
        self.run_code(cmd_name, cmd)

    def hairpin_filter_u2t(self):
        cmd = '{} {} -t '.format(self.python, self.mirna_ref_process)
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
        cmd = '{} {} '.format(self.python, self.mirna_site_filter)
        cmd += '-m {} '.format(self.option('mature_fa').prop['path'])
        cmd += '-s {} '.format(self.mirna_site_json)
        cmd += '-o {}'.format(self.filter_site_json)
        cmd_name = 'mature_filter'
        self.run_code(cmd_name, cmd)

    def mirna_edit_analysis(self):
        cmd = '{} {} '.format(self.python, self.analysis_script)
        cmd += '--hairpin {} '.format(self.hairpin_u2t_hit)
        cmd += '--sample {} '.format(self.filtered_hit)
        cmd += '--site {} '.format(self.filter_site_json)
        cmd += '--name {} '.format(self.mirna_name_json)
        cmd += '--output {} '.format(self.result_xls)
        cmd_name = 'mirna_edit_analysis'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, ignore_error=True):
        command = self.add_command(cmd_name, cmd, ignore_error=ignore_error)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code in [1, -9]:
            self.add_state('memory_limit', 'memory is low!')  # add memory limit error by shicaiping @20200914
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
        source = self.result_xls
        link_name = os.path.join(self.output_dir, '{}.{}'.format(self.option('sample'), os.path.basename(source)))
        if os.path.exists(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('result_xls', link_name)
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
            'id': 'mirna_edit_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'small_rna.mirna_edit',
            'instant': False,
            'options': {
                'species': 'hsa',
                'sample': 'GH1',
                'filtered_fq': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/mirna_edit/tmp/GH1_clip_s.fastq.trimmed',
                'index': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean_index',
                'hairpin_fa': '/mnt/ilustre/users/sanger-dev/app/database/mirbase/hairpin.fa',
                'mature_fa': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/mirna_edit/tmp/mature.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
