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
            dict(name="transcripts_old", type="infile", format="ref_rna_v2.gtf"),
            dict(name="transcripts_new", type="infile", format="ref_rna_v2.gtf"),
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

        self._memory = '20G'

    def end(self):
        super(AsprofileAgent, self).end()


class AsprofileTool(Tool):
    def __init__(self, config):
        super(AsprofileTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python',
            'perl': 'program/perl/perls/perl-5.24.0/bin/perl',
            'extract_as': 'bioinfo/ref_rna_v2/ASprofile/extract-as',
            # 'cat': os.path.join(self.config.SOFTWARE_DIR, 'cat')
        }
        self.script = {
            'count_ref': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/ASprofile/count_ref.py'),
            'extract_as': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/ASprofile/extract-as'),
            'summarize_as': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/ref_rna_v2/ASprofile/summarize_as.pl')
        }
        self.file = {
            'as_old': os.path.join(self.output_dir, '{}_old.as'.format(self.option('sample'))),
            'as_new': os.path.join(self.output_dir, '{}_new.as'.format(self.option('sample'))),
            'old_new': os.path.join(self.output_dir, '{}_AS_result.txt'.format(self.option('sample')))
        }

    def run(self):
        super(AsprofileTool, self).run()
        self.run_extract_as_old()
        # self.run_extract_as_new()
        self.run_summarize_as_old()
        # self.run_summarize_as_new()
        self.run_ref_new()
        # self.run_statistics1()
        self.run_statistics2()
        self.run_statistics3()
        self.set_output()
        self.end()

    # def run_count_ref(self):
    #     cmd = '{} {} {}'.format(self.program['python'], self.script['count_ref'], )

    def run_extract_as_old(self):
        cmd = '{} {} {} > {}'.format(self.program['extract_as'], self.option('transcripts_old').path, self.option('hdrs').path, self.file['as_old'])
        cmd_name = 'run_extract_as_old'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_extract_as_new(self):
        cmd = '{} {} {} > {}'.format(self.program['extract_as'], self.option('transcripts_new').path, self.option('hdrs').path, self.file['as_new'])
        cmd_name = 'run_extract_as_new'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_summarize_as_old(self):
        cmd = '{} {}'.format(self.program['perl'], self.script['summarize_as'])
        cmd += ' {} {}'.format(self.option('transcripts_old').path, self.file['as_old'])
        cmd += ' -p {}_old'.format(self.option('sample'))
        cmd_name = 'run_summarize_as_old'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_summarize_as_new(self):
        cmd = '{} {}'.format(self.program['perl'], self.script['summarize_as'])
        cmd += ' {} {}'.format(self.option('transcripts_new').path, self.file['as_new'])
        cmd += ' -p {}_new'.format(self.option('sample'))
        cmd_name = 'run_summarize_as_new'
        runcmd(self, cmd_name, cmd, shell=True)

    def run_ref_new(self):
        source = os.path.join(self.work_dir, '{}_old.as.nr'.format(self.option('sample')))
        link_name = os.path.join(self.output_dir, '{}_AS_result.txt'.format(self.option('sample')))
        # new_path = os.path.join(self.work_dir, '{}_new.as.nr'.format(self.option('sample')))
        # old = pd.read_table(old_path, sep='\t')
        # old['ref'] = 'old'
        # new = pd.read_table(new_path, sep='\t')
        # new['ref'] = 'new'
        # old_new = pd.concat([old, new], axis=0)
        # old_new.to_csv(self.file['old_new'], index=False, sep='\t')
        os.link(source, link_name)

    # def run_statistics1_old(self):
    #     cmd = 'mv {}/{}_old.as.nr {}/{}_old_AS_result.xls'.format(self.work_dir,
    #                                                       self.option('sample'),
    #                                                       self.output_dir,
    #                                                       self.option('sample'))
    #     cmd_name = 'run_statistics1_old'
    #     runcmd(self, cmd_name, cmd, shell=True)

    # def run_statistics1_new(self):
    #     cmd = 'mv {}/{}_new.as.nr {}/{}_new_AS_result.xls'.format(self.work_dir,
    #                                                       self.option('sample'),
    #                                                       self.output_dir,
    #                                                       self.option('sample'))
    #     cmd_name = 'run_statistics1_new'
    #     runcmd(self, cmd_name, cmd, shell=True)

    def run_statistics2(self):
        os.system('cat {}/{}_AS_result.txt |cut -f 2 |sort |uniq -c > {}/{}_AS_statistics.txt'.format(self.output_dir,
                                                                                           self.option('sample'),
                                                                                           self.output_dir,
                                                                                           self.option('sample')))
        # cmd =  '{} {}/{}_AS_result.txt |cut -f 2 |sort |uniq -c > {}/{}_AS_statistics.txt'.format("cat", self.output_dir,
        #                                                                                    self.option('sample'),
        #                                                                                    self.output_dir,
        #                                                                                    self.option('sample'))

        # command = self.add_command("run_statistics2", cmd, ignore_error=True,shell=True)
        # command.run()
        # self.wait()
        # cmd_name = 'run_statistics2'
        # runcmd(self, cmd_name, cmd, shell=True)

    # def run_statistics2_old(self):
    #     cmd =  'cat {}/{}_old_AS_result.xls |cut -f 2 |sort |uniq -c > {}/{}_old_AS_statistics.xls'.format(self.output_dir,
    #                                                                                        self.option('sample'),
    #                                                                                        self.output_dir,
    #                                                                                        self.option('sample'))
    #     cmd_name = 'run_statistics2_old'
    #     runcmd(self, cmd_name, cmd, shell=True)
    #
    # def run_statistics2_new(self):
    #     cmd =  'cat {}/{}_new_AS_result.xls |cut -f 2 |sort |uniq -c > {}/{}_new_AS_statistics.xls'.format(self.output_dir,
    #                                                                                        self.option('sample'),
    #                                                                                        self.output_dir,
    #                                                                                        self.option('sample'))
    #     cmd_name = 'run_statistics2_new'
    #     runcmd(self, cmd_name, cmd, shell=True)

    def run_statistics3(self):
        path1 = os.path.join(self.output_dir, '{}_AS_statistics.txt'.format(self.option('sample')))
        path2 = os.path.join(self.output_dir, '{}_title_AS_statistics.txt'.format(self.option('sample')))
        event_type = ['AE', 'IR_OFF', 'IR_ON', 'MIR_OFF', 'MIR_ON', 'MSKIP_OFF', 'MSKIP_ON', 'SKIP_OFF', 'SKIP_ON', 'TSS',
                      'TTS', 'XAE', 'XIR_OFF', 'XIR_ON', 'XMIR_OFF', 'XMIR_ON', 'XMSKIP_OFF', 'XMSKIP_ON', 'XSKIP_OFF', 'XSKIP_ON']
        with open(path1, 'r') as file1, open(path2, 'w') as file2:
            file2.write('type' + '\t' + 'num' + '\n')
            for line in file1.readlines():
                number, type = line.strip().split(' ')
                if type == 'event_type':
                    continue
                file2.write(type + '\t' + number + '\n')
                if type in event_type:
                    event_type.remove(type)
            for i in event_type:
                file2.write(i + '\t' + str(0) + '\n')


    # def run_statistics3_new(self):
    #     path1 = os.path.join(self.output_dir, '{}_new_AS_statistics.xls'.format(self.option('sample')))
    #     path2 = os.path.join(self.output_dir, '{}_new_title_AS_statistics.xls'.format(self.option('sample')))
    #     with open(path1, 'r') as file1, open(path2, 'w') as file2:
    #         file2.write('type' + '\t' + 'num' + '\n')
    #         for line in file1.readlines():
    #             number, type = line.strip().split(' ')
    #             if type == 'event_type':
    #                 continue
    #             file2.write(type + '\t' + number + '\n')

    def set_output(self):
        # all ready write results to output
        as_result = os.path.join(self.output_dir, '{}_AS_result.txt'.format(self.option('sample')))
        as_statistics = os.path.join(self.output_dir, '{}_title_AS_statistics.txt'.format(self.option('sample')))
        df = pd.read_table(as_statistics, header=0, sep='\t')
        df.T.to_csv(as_statistics,header=False,index=False,sep='\t')
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

