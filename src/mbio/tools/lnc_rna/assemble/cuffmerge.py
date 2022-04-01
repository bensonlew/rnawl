# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class CuffmergeAgent(Agent):
    '''
    last_modify: 2019.02.11
    '''
    def __init__(self, parent):
        super(CuffmergeAgent, self).__init__(parent)
        options = [
            {'name': 'cpu', 'type': 'int', 'default': 2},
            {'name': 'gtf_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'min_isoform_fraction', 'type': 'float', 'default': 0.1},
            {'name': 'merged_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'merged_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
        ]
        self.add_option(options)
        self.step.add_steps('cuffmerge')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.cuffmerge.start()
        self.step.update()

    def stepfinish(self):
        self.step.cuffmerge.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('cpu', self.option('cpu')))
        if self.option('cpu') == '':
            raise OptionError('number of threads must be specified')
        self.logger.debug('{} - {}'.format('gtf_list', self.option('gtf_list').prop['path']))
        if not self.option('gtf_list').is_set:
            raise OptionError('GTF list file must be provided')
        self.logger.debug('{} - {}'.format('ref_gtf', self.option('ref_gtf').prop['path']))
        if not self.option('ref_gtf').is_set:
            raise OptionError('reference annotation GFF must be provided')
        self.logger.debug('{} - {}'.format('ref_fa', self.option('ref_fa').prop['path']))
        if not self.option('ref_fa').is_set:
            raise OptionError('genomic seqs FASTA must be provided')
        self.logger.debug('{} - {}'.format('min_isoform_fraction', self.option('min_isoform_fraction')))
        if self.option('min_isoform_fraction') == '':
            raise OptionError('minimum isoform fraction must be specified')

    def set_resource(self):
        self._cpu = self.option('cpu')
        self._memory = '16G'

    def end(self):
        super(CuffmergeAgent, self).end()

class CuffmergeTool(Tool):
    def __init__(self, config):
        super(CuffmergeTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/cufflinks-2.2.1'))
        self.cuffmerge = 'bioinfo/rna/cufflinks-2.2.1/cuffmerge'
        self.gffread = 'bioinfo/rna/cufflinks-2.2.1/gffread'
        self.merged_asm = os.path.join(self.work_dir, 'merged_asm')
        self.merged_gtf = os.path.join(self.merged_asm, 'merged.gtf')
        self.merged_fa = os.path.join(self.work_dir, 'merged.fa')

    def run(self):
        super(CuffmergeTool, self).run()
        self.run_cuffmerge()
        self.run_gffread()
        self.set_output()
        self.end()

    def run_cuffmerge(self):
        cmd = '{} {} --keep-tmp'.format(self.cuffmerge, self.option('gtf_list').prop['path'])
        cmd += ' -o {}'.format(self.merged_asm)
        cmd += ' -g {}'.format(self.option('ref_gtf').prop['path'])
        cmd += ' -s {}'.format(self.option('ref_fa').prop['path'])
        cmd += ' --min-isoform-fraction {}'.format(self.option('min_isoform_fraction'))
        cmd += ' -p {}'.format(self.option('cpu'))
        cmd_name = 'run_cuffmerge'
        self.run_code(cmd_name, cmd)

    def run_gffread(self):
        cmd = '{} {}'.format(self.gffread, self.merged_gtf)
        cmd += ' -g {}'.format(self.option('ref_fa').prop['path'])
        cmd += ' -w {}'.format(self.merged_fa)
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
        merged_gtf = os.path.join(self.output_dir, 'merged.gtf')
        if os.path.exists(merged_gtf):
            os.remove(merged_gtf)
        os.link(self.merged_gtf, merged_gtf)
        self.logger.info('succeed in linking {} to {}'.format(self.merged_gtf, merged_gtf))
        self.option('merged_gtf', merged_gtf)
        merged_fa = os.path.join(self.output_dir, 'merged.fa')
        if os.path.exists(merged_fa):
            os.remove(merged_fa)
        os.link(self.merged_fa, merged_fa)
        self.logger.info('succeed in linking {} to {}'.format(self.merged_fa, merged_fa))
        self.option('merged_fa', merged_fa)
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
            'id': 'cuffmerge_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.assemble.cuffmerge',
            'instant': False,
            'options': {
                'gtf_list': '/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/cuffmerge/gtf.list',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()