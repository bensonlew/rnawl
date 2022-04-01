# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class StringtieAgent(Agent):
    '''
    last_modify: 2019.06.26
    '''
    def __init__(self, parent):
        super(StringtieAgent, self).__init__(parent)
        options = [
            {'name': 'fr_stranded', 'type': 'string', 'default': 'fr-unstranded'},
            {'name': 'strand_direct', 'type': 'string', 'default': 'none'},
            {'name': 'sample_name', 'type': 'string', 'default': ''},
            {'name': 'sample_bam', 'type': 'infile', 'format': 'lnc_rna.bam'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'min_coverage', 'type': 'int', 'default': 3},
            {'name': 'min_read', 'type': 'int', 'default': 5},
            {'name': 'cpu', 'type': 'int', 'default': 4},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'sample_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
        ]
        self.add_option(options)
        self.step.add_steps('stringtie')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.stringtie.start()
        self.step.update()

    def stepfinish(self):
        self.step.stringtie.finish()
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
        self.logger.debug('{} - {}'.format('min_coverage', self.option('min_coverage')))
        if self.option('min_coverage') == '':
            raise OptionError('minimum junction coverage must be specified')
        self.logger.debug('{} - {}'.format('min_read', self.option('min_read')))
        if self.option('min_read') == '':
            raise OptionError('minimum reads per bp coverage must be specified')
        self.logger.debug('{} - {}'.format('cpu', self.option('cpu')))
        if self.option('cpu') == '':
            raise OptionError('number of threads must be specified')
        self.logger.debug('{} - {}'.format('ref_fa', self.option('ref_fa').prop['path']))
        if not self.option('ref_fa').is_set:
            raise OptionError('genomic seqs FASTA must be provided')
        self.logger.debug('{} - {}'.format('sample_gtf', self.option('sample_gtf')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('cpu')
        self._memory = '20G'

    def end(self):
        super(StringtieAgent, self).end()

class StringtieTool(Tool):
    def __init__(self, config):
        super(StringtieTool, self).__init__(config)
        self.stringtie = 'bioinfo/rna/stringtie-1.3.4d/stringtie'
        self.gffread = 'bioinfo/rna/cufflinks-2.2.1/gffread'
        self.output_gtf = os.path.join(self.work_dir, '{}_out.gtf'.format(self.option('sample_name')))
        self.output_fa = os.path.join(self.work_dir, '{}_out.fa'.format(self.option('sample_name')))

    def run(self):
        super(StringtieTool, self).run()
        self.run_stringtie()
        self.run_gffread()
        self.set_output()
        self.end()

    def run_stringtie(self):
        if self.option('fr_stranded') == 'fr-unstranded':
            cmd = '{} {} '.format(self.stringtie, self.option('sample_bam').prop['path'])
            cmd += '-G {} '.format(self.option('ref_gtf').prop['path'])
            cmd += '-o {} '.format(self.output_gtf)
            cmd += '-j {} '.format(self.option('min_coverage'))
            cmd += '-c {} '.format(self.option('min_read'))
            cmd += '-p {}'.format(self.option('cpu'))
        elif self.option('fr_stranded') == 'fr-stranded':
            if self.option('strand_direct') == 'firststrand':
                cmd = '{} {} --rf '.format(self.stringtie, self.option('sample_bam').prop['path'])
                cmd += '-G {} '.format(self.option('ref_gtf').prop['path'])
                cmd += '-o {} '.format(self.output_gtf)
                cmd += '-j {} '.format(self.option('min_coverage'))
                cmd += '-c {} '.format(self.option('min_read'))
                cmd += '-p {}'.format(self.option('cpu'))
            elif self.option('strand_direct') == 'secondstrand':
                cmd = '{} {} --fr '.format(self.stringtie, self.option('sample_bam').prop['path'])
                cmd += '-G {} '.format(self.option('ref_gtf').prop['path'])
                cmd += '-o {} '.format(self.output_gtf)
                cmd += '-j {} '.format(self.option('min_coverage'))
                cmd += '-c {} '.format(self.option('min_read'))
                cmd += '-p {}'.format(self.option('cpu'))
        cmd_name = 'run_stringtie'
        self.run_code(cmd_name, cmd)

    def run_gffread(self):
        cmd = '{} {} '.format(self.gffread, self.output_gtf)
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
        for source in [self.output_fa, self.output_gtf]:
            link_name = os.path.join(self.output_dir, os.path.basename(source))
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.option('sample_gtf', link_name)
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
            'id': 'stringtie_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.assemble.stringtie',
            'instant': False,
            'options': {
                'sample_name': '23Y_1115_S6_L002',
                'sample_bam': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/clean_data/23Y_1115_S6_L002.stringtie.sort.bam',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()