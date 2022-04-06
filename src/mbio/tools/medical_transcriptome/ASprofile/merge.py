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


class MergeAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(MergeAgent, self).__init__(parent)
        options = [

            dict(name="result_list", type="infile", format="ref_rna_v2.common"),
            dict(name='statistics_list', type='infile', format='ref_rna_v2.common'),
            dict(name='AS_result_merge', type='outfile', format='ref_rna_v2.common'),
            dict(name='AS_statistics_merge', type='outfile', format='ref_rna_v2.common')

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
        super(MergeAgent, self).end()


class MergeTool(Tool):
    def __init__(self, config):
        super(MergeTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'count_ref': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/ASprofile/count_ref.py'),
        }
        self.file = {
            'result': os.path.join(self.output_dir, 'AS_result_merge.txt'),
            'statistics': os.path.join(self.output_dir, 'AS_statistics_merge.txt')
        }

    def run(self):
        super(MergeTool, self).run()
        self.run_merge()
        self.set_output()
        self.end()

    def run_merge(self):
        result_list = list()
        with open(self.option('result_list').path, 'r') as file:
            for line in file.readlines():
                sample, group, result = line.strip().split('\t')
                df = pd.read_table(result, sep='\t')
                df['sample'] = sample
                df['group'] = group
                result_list.append(df)
        merge_result_matrix = pd.concat(result_list, axis=0)
        merge_result_matrix.to_csv(self.file['result'], index=False, sep='\t')
        statistics_list = list()
        with open(self.option('statistics_list').path, 'r') as file1:
            for line in file1.readlines():
                sample, result = line.strip().split('\t')
                df = pd.read_table(result, sep='\t')
                df['sample'] = sample
                statistics_list.append(df)
        merge_statistics_matrix = pd.concat(statistics_list, axis=0)
        merge_statistics_matrix.to_csv(self.file['statistics'], index=False, sep='\t')


    def set_output(self):
        # all ready write results to output
        self.option('AS_result_merge').set_path(self.file['result'])
        self.option('AS_statistics_merge').set_path(self.file['statistics'])




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

