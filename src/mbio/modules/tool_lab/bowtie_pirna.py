# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from biocluster.module import Module
import sys
import os
import glob
import unittest

class BowtiePirnaModule(Module):
    '''
    last_modify: 2019.04.11
    '''
    def __init__(self, work_id):
        super(BowtiePirnaModule, self).__init__(work_id)
        options = [
            {'name': 'database', 'type': 'infile', 'format': 'ref_rna_v2.common_dir'},
        ]
        self.add_option(options)
        self.bowtie_tools = list()

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
        super(BowtiePirnaModule, self).run()
        self.run_tool()

    def run_tool(self):
        self.bowtie_pirna()
        for tool in self.bowtie_tools:
            tool.run()
        self.on_rely(self.bowtie_tools, self.set_output)

    def bowtie_pirna(self):
        dir = self.option('database').path
        for home, dirs, files in os.walk(dir):
            for filename in files:
                fa_path = os.path.join(home, filename)
                bowtie_pirna = self.add_tool('tool_lab.bowtie_pirna')
                options = {
                    'fa': fa_path,
                }
                bowtie_pirna.set_options(options)
                self.bowtie_tools.append(bowtie_pirna)

    def set_output(self):
        self.end()

    def end(self):
        super(BowtiePirnaModule, self).end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'bowtie_pirna_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'module',
            'name': 'tool_lab.bowtie_pirna',
            'instant': False,
            'options': {
                'database': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/piRNA/ssc/ssc_len/'
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
