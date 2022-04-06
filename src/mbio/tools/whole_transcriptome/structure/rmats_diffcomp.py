# -*- coding: utf-8 -*-
# __author__ = 'gudeqing,qinjincheng'

from biocluster.agent import Agent
from mbio.packages.whole_transcriptome.functions import toolfuncdeco, runcmd
from biocluster.tool import Tool
import os
import unittest

class RmatsDiffcompAgent(Agent):
    '''
    last_modify: 2019.06.14
    '''
    def __init__(self, parent):
        super(RmatsDiffcompAgent, self).__init__(parent)
        options = [
            {'name': 'rmats_detail_list', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'delta_psi', 'type': 'float', 'default': None},
            {'name': 'significant_diff', 'type': 'string', 'default': None},
            {'name': 'significant_value', 'type': 'float', 'default': None},
            {'name': 'diffcomp_txt', 'type': 'outfile', 'format': 'whole_transcriptome.common'},
        ]
        self.add_option(options)
        self.step.add_steps('rmats_diffcomp')
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def step_start(self):
        self.step.rmats_diffcomp.start()
        self.step.update()

    def step_finish(self):
        self.step.rmats_diffcomp.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    @toolfuncdeco
    def end(self):
        super(RmatsDiffcompAgent, self).end()

class RmatsDiffcompTool(Tool):
    def __init__(self,config):
        super(RmatsDiffcompTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'rmats_diffcomp': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/structure/rmats_diffcomp.py')
        }
        self.file = {
            'diffcomp': os.path.join(self.output_dir, 'diffcomp.txt')
        }

    @toolfuncdeco
    def run(self):
        super(RmatsDiffcompTool, self).run()
        self.run_rmats_diffcomp()
        self.set_output()
        self.end()

    @toolfuncdeco
    def run_rmats_diffcomp(self):
        cmd = '{} {}'.format(self.program['python'], self.script['rmats_diffcomp'])
        cmd += ' -i {}'.format(self.option('rmats_detail_list').path)
        cmd += ' -p {}'.format(self.option('delta_psi'))
        cmd += ' -m {}'.format(self.option('significant_diff'))
        cmd += ' -c {}'.format(self.option('significant_value'))
        cmd += ' -o {}'.format(self.file['diffcomp'])
        if self.config.DBVersion == 1:
            DBVersion = 1
        else:
            DBVersion = 0
        cmd += ' -v {}'.format(DBVersion)
        cmd_name = 'run_rmats_diffcomp'
        runcmd(self, cmd_name, cmd)

    @toolfuncdeco
    def set_output(self):
        self.option('diffcomp_txt').set_path(self.file['diffcomp'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test_fdr(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_diffcomp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.structure.rmats_diffcomp',
            'instant': False,
            'options': {
                'rmats_detail_list': '/mnt/ilustre/users/sanger-dev/workspace/20190614/RmatsDiffcomp_tsg_33555_8908_5404/Download/table.list',
                'delta_psi': 0.0,
                'significant_diff': 'fdr',
                'significant_value': 0.05
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

    def test_pvalue(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'rmats_diffcomp_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.structure.rmats_diffcomp',
            'instant': False,
            'options': {
                'rmats_detail_list': '/mnt/ilustre/users/sanger-dev/workspace/20191105/RmatsDiffcomp_workflow_2522_6871/Download/table.list',
                'delta_psi': 0.0,
                'significant_diff': 'pvalue',
                'significant_value': 0.05
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test_pvalue')])
    unittest.TextTestRunner(verbosity=2).run(suite)
