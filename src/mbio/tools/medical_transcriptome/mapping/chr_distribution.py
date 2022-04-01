# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.packages.ref_rna_v3.functions import toolfuncdeco, runcmd


class ChrDistributionAgent(Agent):
    '''
    last_modify: 2019.07.05
    '''

    def __init__(self, parent):
        super(ChrDistributionAgent, self).__init__(parent)
        options = [
            {'name': 'bam', 'type': 'infile', 'format': 'align.bwa.bam'},  # bam格式文件，排序过的
        ]
        self.add_option(options)
        self.step.add_steps('dup')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.dup.start()
        self.step.update()

    def step_end(self):
        self.step.dup.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        if not self.option('bam').is_set:
            raise OptionError('请传入比对结果bam格式文件', code='35601202')
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    @toolfuncdeco
    def end(self):
        super(ChrDistributionAgent, self).end()


class ChrDistributionTool(Tool):
    def __init__(self, config):
        super(ChrDistributionTool, self).__init__(config)
        self.sample = os.path.basename(self.option('bam').path)[:-4]
        self.program = {
            'samtools': 'bioinfo/rna/miniconda2/bin/samtools'
        }
        self.file = {
            'tmp': os.path.join(self.work_dir, '{}.bam_chr_stat.xls'.format(self.sample)),
            'out': os.path.join(self.output_dir, '{}.bam_chr_stat.xls'.format(self.sample))
        }

    @toolfuncdeco
    def run(self):
        super(ChrDistributionTool, self).run()
        self.run_bam2tmp()
        self.run_tmp2out()
        self.set_output()

    @toolfuncdeco
    def run_bam2tmp(self):
        cmd = '{} view {} | cut -f3 | uniq -c > {}'.format(
            self.program['samtools'], self.option('bam').path, self.file['tmp']
        )
        cmd_name = 'run_bam2tmp'
        runcmd(self, cmd_name, cmd, shell=True)

    @toolfuncdeco
    def run_tmp2out(self):
        cmd = r'''awk '{if($2!="*"){print$2"\t"$1}}' %s | ''' % self.file['tmp']
        if len(open(self.file['tmp']).readlines()) >= 31:
            cmd += "awk 'NR==1||$2>100000||length($1)<8' | "
        cmd += 'sort -n -k 1 | '
        cmd += r'''awk 'BEGIN{print"#chr\tread_num"};{print$1"\t"$2}' '''
        cmd += '> %s' % self.file['out']
        cmd_name = 'run_tmp2out'
        runcmd(self, cmd_name, cmd, shell=True)

    @toolfuncdeco
    def set_output(self):
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
            'id': 'chr_distribution_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.mapping.chr_distribution',
            'instant': False,
            'options': {
                'bam': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v3/mouse/bam/A2_2.bam'}
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
