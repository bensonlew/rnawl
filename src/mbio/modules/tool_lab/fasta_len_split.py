# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
import sys
import os
import glob
import unittest

class FastaLenSplitModule(Module):
    '''
    last_modify: 2019.04.11
    '''
    def __init__(self, work_id):
        super(FastaLenSplitModule, self).__init__(work_id)
        options = [
            {'name': 'database', 'type': 'infile', 'format': 'ref_rna_v2.common_dir'},
        ]
        self.add_option(options)
        self.split_tools = list()

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(FastaLenSplitModule, self).run()
        self.run_tool()

    def run_tool(self):
        self.fasta_split()
        for tool in self.split_tools:
            tool.run()
        self.on_rely(self.split_tools, self.set_output)

    def fasta_split(self):
        dir = self.option('database').path
        for home, dirs, files in os.walk(dir):
            if os.path.basename(home) in ['ssc', 'aca']:
                continue
            for filename in files:
                if not filename.endswith('fa'):
                    continue
                fa_path = os.path.join(home, filename)
                len_split = self.add_tool('tool_lab.fasta_len_split')
                options = {
                    'fa': fa_path,
                }
                len_split.set_options(options)
                self.split_tools.append(len_split)

    def set_output(self):
        self.end()

    def end(self):
        super(FastaLenSplitModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'fasta_len_split_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'tool_lab.fasta_len_split',
            'instant': False,
            'options': {
                'database': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/piRNA'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    import sys
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
