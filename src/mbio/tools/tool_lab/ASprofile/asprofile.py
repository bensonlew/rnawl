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


class AsprofileAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(AsprofileAgent, self).__init__(parent)
        options = [

            dict(name="transcripts", type="infile", format="ref_rna_v2.gtf"),
            dict(name="hdrs", type="infile", format="ref_rna_v2.common"),
            dict(name='sample', type='string'),
            dict(name='as_result', type='outfile', format='ref_rna_v2.common'),
            dict(name='as_statistics', type='outfile', format='ref_rna_v2.common')

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
        super(AsprofileAgent, self).end()


class AsprofileTool(Tool):
    def __init__(self, config):
        super(AsprofileTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
            'perl': 'miniconda2/bin/perl',
            'extract_as': 'bioinfo/ref_rna_v2/ASprofile/extract-as'
        }
        self.script = {
            'count_ref': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/ASprofile/count_ref.py'),
            'extract_as': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/ASprofile/extract-as'),
            'summarize_as': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/ASprofile/summarize_as.pl')
        }
        self.file = {
            'as': os.path.join(self.output_dir, '{}.as'.format(self.option('sample'))),
        }

    def run(self):
        super(AsprofileTool, self).run()
        self.run_extract_as()
        self.run_summarize_as()
        self.run_statistics1()
        self.run_statistics2()
        self.run_statistics3()
        self.set_output()
        self.end()

    # def run_count_ref(self):
    #     cmd = '{} {} {}'.format(self.program['python'], self.script['count_ref'], )

    def run_extract_as(self):
        cmd = '{} {} {} > {}'.format(self.program['extract_as'], self.option('transcripts').path, self.option('hdrs').path, self.file['as'])
        cmd_name = 'run_extract_as'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_summarize_as(self):
        cmd = '{} {}'.format(self.program['perl'], self.script['summarize_as'])
        cmd += ' {} {}'.format(self.option('transcripts').path, self.file['as'])
        cmd += ' -p {}'.format(self.option('sample'))
        cmd_name = 'run_summarize_as'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_statistics1(self):
        cmd = 'mv {}/{}.as.nr {}/{}_AS_result.xls'.format(self.work_dir,
                                                          self.option('sample'),
                                                          self.output_dir,
                                                          self.option('sample'))
        cmd_name = 'run_statistics1'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_statistics2(self):
        cmd =  'cat {}/{}_AS_result.xls |cut -f 2 |sort |uniq -c > {}/{}_AS_statistics.xls'.format(self.output_dir,
                                                                                           self.option('sample'),
                                                                                           self.output_dir,
                                                                                           self.option('sample'))
        cmd_name = 'run_statistics2'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_statistics3(self):
        path1 = os.path.join(self.output_dir, '{}_AS_statistics.xls'.format(self.option('sample')))
        path2 = os.path.join(self.output_dir, '{}_new_AS_statistics.xls'.format(self.option('sample')))
        with open(path1, 'r') as file1, open(path2, 'w') as file2:
            file2.write('type' + '\t' + 'num' + '\n')
            for line in file1.readlines():
                number, type = line.strip().split(' ')
                if type == 'event_type':
                    continue
                file2.write(type + '\t' + number + '\n')

    def set_output(self):
        # all ready write results to output
        as_result = os.path.join(self.output_dir, '{}_AS_result.xls'.format(self.option('sample')))
        as_statistics = os.path.join(self.output_dir, '{}_new_AS_statistics.xls'.format(self.option('sample')))
        self.option('as_result').set_path(as_result)
        self.option('as_statistics').set_path(as_statistics)




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "ASprofile" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ASprofile.asprofile",
            "instant": True,
            "options": {
                "transcripts": '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/merged.gtf',
                "hdrs": '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/ref.fa.hdrs',
                'sample': 'Merge'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

