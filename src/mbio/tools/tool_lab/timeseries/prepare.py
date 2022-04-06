# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import unittest

class PrepareAgent(Agent):
    '''
    last_modify: 2019.08.06
    '''
    def __init__(self, parent):
        super(PrepareAgent, self).__init__(parent)
        options = [
            {'name': 'matrix', 'type': 'infile', 'format': 'ref_rna_v2.matrix'},
            {'name': 'geneset', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'pheno', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'method', 'type': 'string', 'default': 'SCM'}, # ['SCM', 'K']
            {'name': 'keep', 'type': 'bool', 'default': False},
            # {'name': 'log', 'type': 'int', 'default': 10},
            # 标准化方法["log_normal"，"normal","no_normal"]
            {"name": "normalized_method", 'type': 'string', 'default': 'log_normal'},
            {'name': 'number', 'type': 'int', 'default': 50},
            {'name': 'unit', 'type': 'int', 'default': 2},
            {'name': 'significance', 'type': 'float', 'default': 0.05},
            {'name': 'correction', 'type': 'string', 'default': 'Bonferroni'}, # ['Bonferroni', 'FDR', 'None']
            {'name': 'clusters', 'type': 'int', 'default': 10},
            {'name': 'starts', 'type': 'int', 'default': 20},
            {'name': 'setting', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps('stem')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.stem.start()
        self.step.update()

    def stepfinish(self):
        self.step.stem.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        self._cpu = 1
        self._memory = '{}G'.format(os.path.getsize(self.option('matrix').path) / 1024 ** 3 + 8)

    @toolfuncdeco
    def end(self):
        super(PrepareAgent, self).end()

class PrepareTool(Tool):
    def __init__(self, config):
        super(PrepareTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'prepare': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/timeseries/prepare.py'),
            'stem': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/stem/stem.jar')
        }
        self.file = {
            'setting': os.path.join(self.output_dir, 'setting.txt')
        }

    @toolfuncdeco
    def run(self):
        super(PrepareTool, self).run()
        self.run_prepare()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_prepare(self):
        cmd = '{} {}'.format(self.program['python'], self.script['prepare'])
        cmd += ' --matrix {}'.format(self.option('matrix').path)
        if self.option('geneset').is_set:
            cmd += ' --geneset {}'.format(self.option('geneset').path)
        cmd += ' --pheno {}'.format(self.option('pheno').path)
        # if self.option('keep'):
        #     cmd += ' --log {}'.format(self.option('log'))
        cmd += ' --normalized_method {}'.format(self.option('normalized_method'))
        cmd += ' --method {}'.format(self.option('method'))
        if self.option('method') == 'SCM':
            cmd += ' --number {}'.format(self.option('number'))
            cmd += ' --unit {}'.format(self.option('unit'))
            cmd += ' --significance {}'.format(self.option('significance'))
            cmd += ' --correction {}'.format(self.option('correction'))
        elif self.option('method') == 'K':
            cmd += ' --clusters {}'.format(self.option('clusters'))
            cmd += ' --starts {}'.format(self.option('starts'))
        cmd += ' --output {}'.format(self.file['setting'])
        cmd_name = 'run_prepare'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def set_output(self):
        self.option('setting').set_path(self.file['setting'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_scm(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'prepare_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.timeseries.prepare',
            'instant': False,
            'options': {
                'matrix': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/matrix.tsv',
                'geneset': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/geneset.list',
                'pheno': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/ref_rna_v2/stem/pheno.txt',
                'method': 'SCM',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_scm')])
    unittest.TextTestRunner(verbosity=2).run(suite)
