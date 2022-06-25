# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.packages.ref_rna_v3.functions import toolfuncdeco, runcmd


class CoverageAgent(Agent):
    """
    last_modify: 2019.11.07
    """

    def __init__(self, parent):
        super(CoverageAgent, self).__init__(parent)
        options = [
            {'name': 'bam', 'type': 'infile', 'format': 'align.bwa.bam'},  # bam格式文件，排序过的
            {'name': 'bed', 'type': 'infile', 'format': 'ref_rna_v2.bed'},  # bed格式文件
            {'name': 'min_len', 'type': 'int', 'default': 100}  # Minimum mRNA length (bp).
        ]
        self.add_option(options)
        self._memory_increase_step = 20
        self.step.add_steps('coverage')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.coverage.start()
        self.step.update()

    def step_end(self):
        self.step.coverage.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        if not self.option('bam').is_set:
            raise OptionError('请传入比对结果bam格式文件', code='35600303')
        if not self.option('bed').is_set:
            raise OptionError('请传入bed格式文件', code='35600304')
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    @toolfuncdeco
    def end(self):
        super(CoverageAgent, self).end()


class CoverageTool(Tool):
    def __init__(self, config):
        super(CoverageTool, self).__init__(config)
        self.prefix = os.path.join(self.work_dir, os.path.basename(self.option('bam').path))[:-4]
        self.program = {
            'samtools': 'miniconda2/bin/samtools',
            'python': 'miniconda2/bin/python',
            'bash': 'sh'
        }
        self.script = {
            'genebody_coverage': os.path.join(self.config.SOFTWARE_DIR,
                                              'miniconda2/bin/geneBody_coverage.py')
        }
        self.file = {
            'bam': '{}.bam'.format(self.prefix),
            'sh': os.path.join(self.work_dir, 'coverage.sh'),
            'txt': '{}.geneBodyCoverage.txt'.format(self.prefix)
        }

    def run(self):
        super(CoverageTool, self).run()
        if not self.check_completed():
            self.run_samtools_index()
            self.run_bash()
        self.set_output()

    def check_completed(self):
        if os.path.isfile(self.file['txt']):
            lines = open(self.file['txt']).readlines()
        else:
            return False
        if len(lines) == 2:
            for line in lines:
                if len(line.strip().split('\t')) == 101:
                    continue
                else:
                    break
            else:
                return True
        else:
            return False

    def run_samtools_index(self):
        if os.path.isfile(self.file['bam']):
            os.remove(self.file['bam'])
        os.link(self.option('bam').path, self.file['bam'])
        cmd = '{} index {}'.format(self.program['samtools'], self.file['bam'])
        runcmd(self, 'run_samtools_index', cmd)

    def run_genebody_coverage(self):
        cmd = '{} {}'.format(self.program['python'], self.script['genebody_coverage'])
        cmd += ' -i {}'.format(self.file['bam'])
        cmd += ' -r {}'.format(self.option('bed').path)
        cmd += ' -l {}'.format(self.option('min_len'))
        cmd += ' -o {}'.format(self.prefix)
        runcmd(self, 'run_genebody_coverage', cmd)

    def run_bash(self):
        text = '{} {}'.format(os.path.join(self.config.SOFTWARE_DIR, self.program['python']),
                              self.script['genebody_coverage'])
        text += ' -i {}'.format(self.file['bam'])
        text += ' -r {}'.format(self.option('bed').path)
        text += ' -l {}'.format(self.option('min_len'))
        text += ' -o {}'.format(self.prefix)
        open(self.file['sh'], 'w').write(text)
        cmd = '{} {}'.format(self.program['bash'], self.file['sh'])
        runcmd(self, 'run_bash', cmd, shell=True)

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
            'id': 'coverage_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.mapping.coverage',
            'instant': False,
            'options': {
                'bam': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v3/coverage/B3_1.bam',
                'bed': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v3/coverage/Brapa_genome_v3.0.gtf.filter.bed',
                'min_len': 100
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
