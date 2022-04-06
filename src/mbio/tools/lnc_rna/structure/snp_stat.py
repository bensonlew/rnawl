# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong, qinjincheng'

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
import os
import unittest

class SnpStatAgent(Agent):
    '''
    last_modify: 2019.03.06
    '''
    def __init__(self, parent):
        super(SnpStatAgent, self).__init__(parent)
        options = [
            {'name': 'snp_detail', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'snp_annotation', 'type': 'infile', 'format': 'lnc_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps('snp_stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.snp_stat.start()
        self.step.update()

    def step_end(self):
        self.step.snp_stat.finish()
        self.step.update()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        if self.option('snp_detail').is_set:
            self.logger.debug('{} = {}'.format('snp_detail', self.option('snp_detail').path))
            self.infile_size = os.path.getsize(self.option('snp_detail').path)
        else:
            raise OptionError('SNP detail file must be provided')
        if self.option('snp_annotation').is_set:
            self.logger.debug('{} = {}'.format('snp_annotation', self.option('snp_annotation').path))
        else:
            raise OptionError('SNP annotation file must be provided')

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(int(float(self.infile_size) / 1024 ** 3 * 48 + 16))

    def end(self):
        super(SnpStatAgent, self).end()

class SnpStatTool(Tool):
    def __init__(self, config):
        super(SnpStatTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.snp_statistics_py = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/snp_statistics.py')

    def run(self):
        super(SnpStatTool, self).run()
        self.run_snp_statistics()
        self.set_output()
        self.end()

    def run_snp_statistics(self):
        cmd = '{} {}'.format(self.python, self.snp_statistics_py)
        cmd += ' -d {}'.format(self.option('snp_detail').path)
        cmd += ' -a {}'.format(self.option('snp_annotation').path)
        cmd += ' -o {}'.format(self.work_dir)
        cmd_name = 'run_snp_statistics'
        self.run_code(cmd_name, cmd)

    def run_code(self, cmd_name, cmd, shell=False):
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
        basenames = [
            'data_anno_pre.xls',
            'indel_position_distribution.xls',
            'main_info.txt',
            'snp_annotation_statistics.xls',
            'snp_depth_statistics.xls',
            'snp_freq_statistics.xls',
            'snp_position_distribution.xls',
            'snp_transition_tranversion_statistics.xls',
        ]
        for basename in basenames:
            source = os.path.join(self.work_dir, basename)
            link_name = os.path.join(self.output_dir, basename)
            if os.path.exists(link_name):
                os.remove(link_name)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
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
            'id': 'snp_stat_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'lnc_rna.structure.snp_stat',
            'instant': False,
            'options': {
                'snp_detail': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/annovar/output/snp_detail.xls',
                'snp_annotation': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/annovar/output/snp_annotation.xls',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
