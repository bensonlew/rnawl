# -*- coding: utf-8 -*-
# __author__ = 'zengjing,qinjincheng'

import glob
import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.packages.ref_rna_v3.functions import toolfuncdeco, runcmd


class BamReadsdistributionAgent(Agent):
    '''
    last_modify: 2019.07.05
    '''

    def __init__(self, parent):
        super(BamReadsdistributionAgent, self).__init__(parent)
        options = [
            {'name': 'bam', 'type': 'infile', 'format': 'align.bwa.bam'},  # bam格式文件，排序过的
            {'name': 'bed', 'type': 'infile', 'format': 'ref_rna_v2.bed'}  # bed格式文件
        ]
        self.add_option(options)
        self._memory_increase_step = 30
        self.step.add_steps('reads_distribution')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.reads_distribution.start()
        self.step.update()

    def step_end(self):
        self.step.reads_distribution.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        if not self.option('bam').is_set:
            raise OptionError('请传入比对结果bam格式文件', code='35600503')
        if not self.option('bed').is_set:
            raise OptionError('请传入参考基因组结构注释bed格式文件', code='35600504')
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 2
        self._memory = '40G'

    @toolfuncdeco
    def end(self):
        super(BamReadsdistributionAgent, self).end()


class BamReadsdistributionTool(Tool):
    def __init__(self, config):
        super(BamReadsdistributionTool, self).__init__(config)
        self.prefix = os.path.join(self.work_dir, os.path.basename(self.option('bam').path)[:-4])
        self.program = {
            'python': self.config.SOFTWARE_DIR + '/miniconda2/bin/python'
        }
        self.script = {
            'read_distribution': os.path.join(self.config.SOFTWARE_DIR,
                                              'miniconda2/bin/read_distribution.py')
        }
        self.file = {
            'txt': '{}.reads_distribution.txt'.format(self.prefix)
        }
        self.cmd_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.shell_path = 'bioinfo/rna/scripts'

    @toolfuncdeco
    def run(self):
        super(BamReadsdistributionTool, self).run()
        self.run_reads_distribution()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_reads_distribution(self):
        cmd = '{} {}'.format(self.program['python'], self.script['read_distribution'])
        cmd += ' -i {}'.format(self.option('bam').path)
        cmd += ' -r {}'.format(self.option('bed').path)
        cmd += ' > {}'.format(self.file['txt'])
        cmd_name = 'run_reads_distribution'
        # runcmd(self, cmd_name, cmd, shell=True)
        command = self.add_command(cmd_name, cmd, shell=True, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code in [137, 139]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")

    @toolfuncdeco
    def set_output(self):
        source = self.file['txt']
        link_name = os.path.join(self.output_dir, os.path.basename(source))
        if os.path.isfile(link_name):
            os.remove(link_name)
        os.link(source, link_name)
        self.end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bam_readsdistribution_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.mapping.bam_readsdistribution',
            'instant': False,
            'options': {
                'bam': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/sanger_test/RT_EB2.bam',
                'bed': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/sanger_test/final.gtf.filter.bed'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
