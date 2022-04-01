# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class PalindromeAgent(Agent):
    '''
    last_modify: 2020.07.23
    '''

    def __init__(self, parent):
        super(PalindromeAgent, self).__init__(parent)
        options = [
            {'name': 'fa_input', 'type': 'infile', 'format': 'whole_transcriptome.fasta'},
            {'name': 'range_min', 'type': 'int', 'default': 6},
            {'name': 'range_max', 'type': 'int', 'default': 30},
            {'name': 'palindrome_table', 'type': 'outfile', 'format': 'whole_transcriptome.common'}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(PalindromeAgent, self).end()


class PalindromeTool(Tool):
    def __init__(self, config):
        super(PalindromeTool, self).__init__(config)
        self.program = {
            'python': os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin/python')
        }
        self.script = {
            'palindrome': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/palindrome.py')
        }
        self.file = {
            'palindrome_table': os.path.join(self.output_dir, "palindrome.txt")
        }

    def run(self):
        super(PalindromeTool, self).run()
        self.run_palindrome()
        self.set_output()
        self.end()

    def run_palindrome(self):
        cmd = '{} {} -f {} -mi {} -ma {} -o {}'.format(self.program['python'], self.script['palindrome'],
                                                       self.option('fa_input').path, self.option('range_min'),
                                                       self.option('range_max'), self.file['palindrome_table'])
        cmd_name = 'run_palindrom'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
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

    def set_output(self):
        self.option('palindrome_table').set_path(self.file['palindrome_table'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'palindrome_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.palindrome',
            'instant': False,
            'options': {
                'fa_input': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/palindrome/test.fasta',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


