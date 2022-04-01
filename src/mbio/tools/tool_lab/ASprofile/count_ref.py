# coding=utf-8
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
__author__ = 'zoujiaxun'


class CountRefAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(CountRefAgent, self).__init__(parent)
        options = [

            dict(name="ref_fa", type="infile", format="ref_rna_v2.fasta"),
            dict(name='hdrs', type='outfile', format='ref_rna_v2.common')


        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        pass

    def set_resource(self):
        self._cpu = 1

        self._memory = '10G'

    def end(self):
        super(CountRefAgent, self).end()


class CountRefTool(Tool):
    def __init__(self, config):
        super(CountRefTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
        }
        self.script = {
            'count_ref': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/ASprofile/count_ref.py'),
        }
        self.file = {
            'hdrs': os.path.join(self.output_dir, 'ref.fa.hdrs'),
        }

    def run(self):
        super(CountRefTool, self).run()
        self.run_count_ref()
        self.set_output()
        self.end()

    def run_count_ref(self):
        cmd = '{} {} {} {}'.format(self.program['python'], self.script['count_ref'], self.option('ref_fa').path, self.file['hdrs'])
        cmd_name = 'run_count_ref'
        runcmd(self, cmd_name, cmd)


    def set_output(self):
        # all ready write results to output
        self.option('hdrs').set_path(self.file['hdrs'])




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "count_ref" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ASprofile.count_ref",
            "instant": True,
            "options": {
                "ref_fa": '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/ref.fa',

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

