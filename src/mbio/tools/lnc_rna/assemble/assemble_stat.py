# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue, shicaiping, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class AssembleStatAgent(Agent):
    '''
    last_modify: 2019.01.27
    '''
    def __init__(self, parent):
        super(AssembleStatAgent, self).__init__(parent)
        options = [
            {'name': 'all_files_dir', 'type': 'infile', 'format': 'lnc_rna.common_dir'},
            {'name': 'assemble_method', 'type': 'string', 'default': ''},
        ]
        self.add_option(options)
        self.step.add_steps('assemble_stat')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.assemble_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.assemble_stat.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        self.logger.debug('{} - {}'.format('all_files_dir', self.option('all_files_dir').prop['path']))
        if not self.option('all_files_dir').is_set:
            raise OptionError('directory containing all result files must be provided')
        self.logger.debug('{} - {}'.format('assemble_method', self.option('assemble_method')))
        if self.option('assemble_method') == '':
            raise OptionError('assemble_method must be specified')
        self.logger.info('finish check_options at {}'.format(self.__class__.__name__))

    def set_resource(self):
        self._cpu = 1
        self._memory = '8G'

    def end(self):
        super(AssembleStatAgent, self).end()

class AssembleStatTool(Tool):
    def __init__(self, config):
        super(AssembleStatTool, self).__init__(config)
        self.python = 'program/Python/bin/python'
        self.get_number_list_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/get_number_list.py')
        self.step_count_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/step_count.py')
        self.class_code_count_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/class_code_count.py')
        self.gene_trans_exon_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/gene_trans_exon.py')
        self.count_trans_exon_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/count_trans_exon.py')
        self.number_list = os.path.join(self.work_dir, 'number_list.txt')
        self.all_transcripts_fa = os.path.join(self.option('all_files_dir').prop['path'], 'all_transcripts.fa')
        self.merged_fa_txt = os.path.join(self.option('all_files_dir').prop['path'], 'merged.fa.txt')
        self.add_code_merged_gtf = os.path.join(self.option('all_files_dir').prop['path'], 'add_code_merged.gtf')
        self.code_num_txt = os.path.join(self.option('all_files_dir').prop['path'], 'code_num.txt')

    def run(self):
        super(AssembleStatTool, self).run()
        self.run_get_number_list()
        self.run_step_count()
        self.run_class_code_count()
        self.run_gene_trans_exon()
        self.run_count_trans_exon()
        self.set_output()
        self.end()

    def run_get_number_list(self):
        '''
        import *_out.gtf and merged.gtf
        export number_list.txt
        '''
        cmd = '{} {}'.format(self.python, self.get_number_list_py)
        cmd += ' -i {}'.format(self.option('all_files_dir').prop['path'])
        cmd += ' -o {}'.format(self.number_list)
        cmd_name = 'run_get_number_list'
        self.run_code(cmd_name, cmd)

    def run_step_count(self):
        '''
        import merged.fa
        export trans_count_stat_*.txt
        '''
        cmd = '{} {}'.format(self.python, self.step_count_py)
        cmd += ' -f {}'.format(self.all_transcripts_fa)
        cmd += ' -d {}'.format(self.option('all_files_dir').prop['path'])
        cmd += ' -g {}'.format('10')
        cmd += ' -s {}'.format('200,300,600,1000')
        cmd += ' -o {}'.format(self.work_dir)
        cmd_name = 'run_step_count'
        self.run_code(cmd_name, cmd)

    def run_class_code_count(self):
        '''
        import add_code_merged.gtf
        export code_num.txt
        '''
        cmd = '{} {}'.format(self.python, self.class_code_count_py)
        cmd += ' -i {}'.format(self.add_code_merged_gtf)
        cmd += ' -o {}'.format(self.code_num_txt)
        cmd_name = 'run_class_code_count'
        self.run_code(cmd_name, cmd)

    def run_gene_trans_exon(self):
        '''
        import old_*.gtf and new_*.gtf
        export old_*.gtf.trans, old_*.gtf.exon, new_*.gtf.trans, new_*.gtf.exon
        '''
        cmd = '{} {}'.format(self.python, self.gene_trans_exon_py)
        cmd += ' -i {}'.format(self.option('all_files_dir').prop['path'])
        cmd += ' -m {}'.format(self.option('assemble_method'))
        cmd += ' -o {}'.format(self.work_dir)
        cmd_name = 'run_gene_trans_exon'
        self.run_code(cmd_name, cmd)

    def run_count_trans_exon(self):
        '''
        import *.gtf.trans and *.gtf.exon
        export *.gtf.trans_*.txt and *.gtf.exon_*.txt
        '''
        cmd = '{} {}'.format(self.python, self.count_trans_exon_py)
        cmd += ' -i {}'.format(self.work_dir)
        cmd += ' -s {}'.format('1,5,10,20')
        cmd += ' -o {}'.format(self.work_dir)
        cmd_name = 'run_count_trans_exon'
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
        for f in os.listdir(self.work_dir):
            if f.endswith('.txt') and f != 'log.txt' and not f.endswith('_resource.txt'):
                source = os.path.join(self.work_dir, f)
                link_name = os.path.join(self.output_dir, f)
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
                self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        code_num_txt = os.path.join(self.output_dir, os.path.basename(self.code_num_txt))
        if os.path.exists(code_num_txt):
            os.remove(code_num_txt)
        os.link(self.code_num_txt, code_num_txt)
        self.logger.info('succeed in linking {} to {}'.format(self.code_num_txt, code_num_txt))
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
            'id': 'assemble_stat_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.assemble.assemble_stat',
            'instant': False,
            'options': {
                'all_files_dir': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/assemble',
                'assemble_method': 'stringtie',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()