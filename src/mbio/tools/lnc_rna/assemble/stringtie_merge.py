# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class StringtieMergeAgent(Agent):
    '''
    last_modify: 2019.01.24
    '''
    def __init__(self, parent):
        super(StringtieMergeAgent, self).__init__(parent)
        options = [
            {'name': 'cpu', 'type': 'int', 'default': 2},
            {'name': 'gtf_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'min_cov', 'type': 'int', 'default': 5},
            {'name': 'min_tpm', 'type': 'int', 'default': 1},
            {'name': 'min_iso', 'type': 'float', 'default': 0.1},
            {'name': 'ref_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'merged_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'merged_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
        ]
        self.add_option(options)
        self.step.add_steps('stringtie_merge')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.stringtie_merge.start()
        self.step.update()

    def stepfinish(self):
        self.step.stringtie_merge.finish()
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
        self.logger.debug('{} - {}'.format('min_cov', self.option('min_cov')))
        if self.option('min_cov') == '':
            raise OptionError('minimum input transcript coverage must be specified')
        self.logger.debug('{} - {}'.format('min_tpm', self.option('min_tpm')))
        if self.option('min_tpm') == '':
            raise OptionError('minimum input transcript TPM must be specified')
        self.logger.debug('{} - {}'.format('min_iso', self.option('min_iso')))
        if self.option('min_iso') == '':
            raise OptionError('minimum isoform fraction must be specified')
        self.logger.debug('{} - {}'.format('ref_fa', self.option('ref_fa').prop['path']))
        if not self.option('ref_fa').is_set:
            raise OptionError('genomic seqs FASTA must be provided')
        self.logger.debug('{} - {}'.format('merged_gtf', self.option('merged_gtf')))
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = self.option('cpu')
        self._memory = '8G'

    def end(self):
        super(StringtieMergeAgent, self).end()

class StringtieMergeTool(Tool):
    def __init__(self, config):
        super(StringtieMergeTool, self).__init__(config)
        self.stringtie = 'bioinfo/rna/stringtie-1.3.4d/stringtie'
        self.gffread = 'bioinfo/rna/cufflinks-2.2.1/gffread'
        self.python = 'program/Python/bin/python'
        self.filter_gtf_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/filter_gtf.py')
        self.merge_raw_gtf = os.path.join(self.work_dir, 'merge.raw.gtf')
        self.merge_filter_gtf = os.path.join(self.work_dir, 'merge.filter.gtf')
        self.merge_filter_fa = os.path.join(self.work_dir, 'merge.filter.fa')

    def run(self):
        super(StringtieMergeTool, self).run()
        self.reset_gtflist()
        self.run_stringtie_merge()
        self.run_filter_gtf()
        self.run_gffread()
        self.set_output()
        self.end()

    def run_stringtie_merge(self):
        cmd = '{} --merge {} '.format(self.stringtie, 'gtf_clean.list')
        cmd += '-G {} '.format(self.option('ref_gtf').prop['path'])
        cmd += '-o {} '.format(self.merge_raw_gtf)
        cmd += '-c {} '.format(self.option('min_cov'))
        cmd += '-T {} '.format(self.option('min_tpm'))
        cmd += '-f {} '.format(self.option('min_iso'))
        cmd_name = 'run_stringtie_merge'
        self.run_code(cmd_name, cmd)

    def run_filter_gtf(self):
        cmd = '{} {} '.format(self.python, self.filter_gtf_py)
        cmd += '--ref {} '.format(self.option('ref_gtf').prop['path'])
        cmd += '--raw {} '.format(self.merge_raw_gtf)
        cmd += '--filter {}'.format(self.merge_filter_gtf)
        cmd_name = 'run_filter_gtf'
        self.run_code(cmd_name, cmd)

    def run_gffread(self):
        cmd = '{} {} '.format(self.gffread, self.merge_filter_gtf)
        cmd += '-g {} '.format(self.option('ref_fa').prop['path'])
        cmd += '-w {}'.format(self.merge_filter_fa)
        cmd_name = 'run_gffread'
        self.run_code(cmd_name, cmd)

    def reset_gtflist(self):
        '''
        删除空gtf
        '''
        with open(self.option('gtf_list').prop['path'], 'r') as f, open("gtf_clean.list", 'w') as fo:
            for line in f:
                with open(line.strip(), 'r') as f:
                    if len(f.readlines()) >= 3:
                        fo.write(line)

    def run_code(self, cmd_name, cmd):
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('succeed in running {}'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('fail to run {}, try again'.format(cmd_name))
            '''
            if cmd_name == "run_stringtie_merge" and command.return_code is 1:
                self.logger.warn('gtf 中存在空文件?')
            
            self.reset_gtflist()
            '''
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
        os.link(self.merge_filter_gtf, merged_gtf)
        self.logger.info('succeed in linking {} to {}'.format(self.merge_filter_gtf, merged_gtf))
        self.option('merged_gtf', merged_gtf)
        merged_fa = os.path.join(self.output_dir, 'merged.fa')
        if os.path.exists(merged_fa):
            os.remove(merged_fa)
        os.link(self.merge_filter_fa, merged_fa)
        self.logger.info('succeed in linking {} to {}'.format(self.merge_filter_gtf, merged_fa))
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
            'id': 'stringtie_merge_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.assemble.stringtie_merge',
            'instant': False,
            'options': {
                'gtf_list': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/assemble/stringtie.gtf.list',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
