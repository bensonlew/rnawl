# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.ref_rna_v2.functions import runcmd


class MasigproAgent(Agent):
    '''
    last_modify: 2019.12.18
    '''

    def __init__(self, parent):
        super(MasigproAgent, self).__init__(parent)
        options = [
            {'name': 'matrix', 'type': 'infile', 'format': 'ref_rna_v2.matrix'},
            {'name': 'geneset', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'design', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'cluster', 'type': 'int', 'default': 9},
            {'name': 'method', 'type': 'string', 'default': 'hclust'},  # ['hclust', 'kmeans', 'Mclust']
            {'name': 'level', 'type': 'string', 'default': 'T'},  # ['T', 'G']
            {'name': 'way', 'type': 'string', 'default': 'tpm'},  # ['tpm', 'fpkm']
            {'name': 'result', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps('masigpro')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.masigpro.start()
        self.step.update()

    def stepfinish(self):
        self.step.masigpro.finish()
        self.step.update()

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('matrix').path) / 1024 ** 3 + 8)

    def end(self):
        super(MasigproAgent, self).end()


class MasigproTool(Tool):
    def __init__(self, config):
        super(MasigproTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': 'bioinfo/rna/miniconda2/bin/Rscript'
        }
        self.script = {
            'pretreat': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/timeseries/pretreat.py'),
            'masigpro': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v2/timeseries/masigpro.r'),
            'heatmap.py': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/timeseries/heatmap.py'),
            'heatmap.r': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/timeseries/heatmap.r')
        }
        self.file = {
            'matrix': os.path.join(self.work_dir, 'matrix.delim.txt'),
            'desgin': os.path.join(self.work_dir, 'design.delim.txt'),
            'result': os.path.join(self.output_dir, 'result.tsv')
        }

    def run(self):
        super(MasigproTool, self).run()
        self.run_pretreat()
        self.run_masigpro()
        self.run_heatmap()
        self.set_output()
        self.end()

    def run_pretreat(self):
        cmd = '{} {}'.format(self.program['python'], self.script['pretreat'])
        cmd += ' -m {}'.format(self.option('matrix').path)
        if self.option('geneset').is_set:
            cmd += ' -g {}'.format(self.option('geneset').path)
        cmd += ' -d {}'.format(self.option('design').path)
        cmd += ' -l {}'.format(self.option('level'))
        cmd += ' -o {}'.format(self.work_dir)
        cmd_name = 'run_pretreat'
        runcmd(self, cmd_name, cmd)

    def run_masigpro(self):
        cmd = '{} {}'.format(self.program['rscript'], self.script['masigpro'])
        cmd += ' -i {}'.format(self.file['matrix'])
        cmd += ' -d {}'.format(self.file['desgin'])
        cmd += ' -k {}'.format(self.option('cluster'))
        cmd += ' -m {}'.format(self.option('method'))
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_masigpro'
        runcmd(self, cmd_name, cmd)

    def run_heatmap(self):
        cmd = '{} {}'.format(self.program['python'], self.script['heatmap.py'])
        cmd += ' -m {}'.format(self.option('matrix').path)
        cmd += ' -d {}'.format(self.option('design').path)
        cmd += ' -r {}'.format(self.file['result'])
        cmd += ' --interpreter {}'.format(os.path.join(self.config.SOFTWARE_DIR, self.program['rscript']))
        cmd += ' --script {}'.format(self.script['heatmap.r'])
        cmd += ' -l {}'.format(self.option('level'))
        cmd += ' -t {}'.format(self.option('way'))
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_heatmap'
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        self.option('result').set_path(self.file['result'])


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test_hclust(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'masigpro_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.timeseries.masigpro',
            'instant': False,
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/matrix.tsv',
                'geneset': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/geneset.list',
                'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/design.txt',
                'cluster': 9,
                'method': 'hclust'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_kmeans(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'masigpro_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.timeseries.masigpro',
            'instant': False,
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.matrix.tsv',
                'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.design.tsv',
                'cluster': 4,
                'method': 'kmeans'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_mclust(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'masigpro_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.timeseries.masigpro',
            'instant': False,
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.matrix.tsv',
                'design': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/masigpro/vs.design.tsv',
                'cluster': 4,
                'method': 'Mclust'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_kmeans')])
    unittest.TextTestRunner(verbosity=2).run(suite)
