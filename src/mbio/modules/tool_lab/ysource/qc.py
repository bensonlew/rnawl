# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import unittest
import random
import datetime
import glob
import subprocess
import re
import os
import sys
import shutil


class QcModule(Module):
    def __init__(self, work_id):
        super(QcModule, self).__init__(work_id)
        options = [
            {"name": "sample_list", "type": "string"},
            {"name": "quality_score", "type": "int", "default": 20},
            {"name": "min_length", "type": "int", "default": 25},
            {"name": "f_adapter", "type": "string",
                "default": "GATCGGAAGAGCACACGTCT"},
            {"name": "r_adapter", "type": "string",
                "default": "AATGATACGGCGACCACCGA"},
            # {"name":"quality_score","type":"int","default":20},
            # {"name":"min_length","type":"int","default": 25},
            {"name": "trim_ns", "type": "int", "default": 1},
            {"name": "trim_left", "type": "int", "default": 0},
            {"name": "ns_max_n", "type": "int", "default": 5},
            {"name": "trim_qual_window", "type": "int", "default": 5},
            {"name": "trim_qual_step", "type": "int", "default": 1},
            {"name": "line_width", "type": "int", "default": 1}
        ]
        self.add_option(options)
        self.sp_info = {}
        self.sp_list = []
        self.prinseq_tools = list()
        self.seqprep_tools = list()
        self.prinseq_option = {}
        # self.prinseq = self.add_tool('tool_lab.ysource.prinseq')
        # self.seqprep = self.add_tool('tool_lab.ysource.seqprep')

    def check_options(self):
        if not os.path.exists(self.option("sample_list")):
            raise OptionError("请输入sample_list文件")
        with open(self.option('sample_list'), 'r') as snl:
            while 1:
                line = snl.readline()
                self.logger.info(line)
                if not line:
                    break
                field = line.rstrip().split('\t')
                sn = field[1]
                fq1 = field[2]
                fq2 = field[3]
                self.sp_list.append(sn)
                if not os.path.exists(fq1) and not os.path.exists(fq2):
                    raise OptionError("缺少fq.gz文件")
                self.sp_info[sn] = [fq1, fq2]
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(QcModule, self).run()
        self.run_tool()

    def run_tool(self):
        self.run_seqprep()
        self.logger.info(len(self.seqprep_tools))
        if len(self.seqprep_tools) == 1:
            self.seqprep_tools[0].on('end', self.run_prinseq)
        else:
            self.on_rely(self.seqprep_tools, self.run_prinseq)
        for tool in self.seqprep_tools:
            tool.run()

    def run_seqprep(self):
        for i in self.sp_list:
            sample_name = i
            fq1 = self.sp_info[i][0]
            fq2 = self.sp_info[i][1]
            seqprep = self.add_tool('tool_lab.ysource.seqprep')
            options = {
                "sample_name": sample_name,
                "R1": fq1,
                "R2": fq2,
                "quality_score": self.option('quality_score'),
                "min_length": self.option("min_length"),
                "f_adapter": self.option("f_adapter"),
                "r_adapter": self.option("r_adapter"),
            }
            self.prinseq_option[i] = {
                "sample_name": i,
                "f1": seqprep.output_dir + "/R1_cutadapt_{}".format(i),
                "f2": seqprep.output_dir + "/R2_cutadapt_{}".format(i),
                "quality_score": self.option('quality_score'),
                "min_length": self.option("min_length"),
                "trim_ns": self.option("trim_ns"),
                "trim_left": self.option("trim_left"),
                "ns_max_n": self.option("ns_max_n"),
                "trim_qual_window": self.option("trim_qual_window"),
                "trim_qual_step": self.option("trim_qual_step"),
                "line_width": self.option("line_width")
            }
            seqprep.set_options(options)
            self.seqprep_tools.append(seqprep)

    def run_prinseq(self):
        for i in self.sp_list:
            sample_name = i
            index = self.sp_list.index(i)
            prinseq = self.add_tool('tool_lab.ysource.prinseq')
            # options = {
            #     "sample_name": i,
            #     "f1": self.seqprep_tools[index].output_dir + "/R1_cutadapt_{}".format(i),
            #     "f2": self.seqprep_tools[index].output_dir + "/R2_cutadapt_{}".format(i),
            #     "quality_score": self.option('quality_score'),
            #     "min_length": self.option("min_length"),
            #     "trim_ns": self.option("trim_ns"),
            #     "trim_left": self.option("trim_left"),
            #     "ns_max_n": self.option("ns_max_n"),
            #     "trim_qual_window": self.option("trim_qual_window"),
            #     "trim_qual_step": self.option("trim_qual_step"),
            #     "line_width": self.option("line_width")
            # }
            prinseq.set_options(self.prinseq_option[i])
            self.prinseq_tools.append(prinseq)
        if len(self.prinseq_tools) == 1:
            self.prinseq_tools[0].on('end', self.set_output)
        else:
            self.on_rely(self.prinseq_tools, self.set_output)
        for tool in self.prinseq_tools:
            tool.run()

    def set_output(self, event):
        self.logger.info('start set_output at {}'.format(
            self.__class__.__name__))
        with open(os.path.join(self.output_dir, "bam_list"), 'w') as bl:
            for tool in self.prinseq_tools:
                sample_name = ""
                text1 = ""
                text2 = ""
                for source in glob.glob(os.path.join(tool.output_dir, '*')):
                    clean_input = open(os.path.join(self.output_dir,"{}_clean_input.txt".format(tool.option('sample_name'))),"w")
                    link_name = os.path.join(
                        self.output_dir, os.path.basename(source))
                    if os.path.exists(link_name):
                        os.remove(link_name)
                    os.link(source, link_name)
                    sample_name = tool.option('sample_name')
                    if os.path.basename(source)[-7] == "1":
                        text1 = link_name
                    elif os.path.basename(source)[-7] == "2":
                        text2 = link_name
                    clean_input.write("CLEAN\t"+ text1 +'\t'+ text2)
                    clean_input.close()
                    # text1 = text1 + "\t" + link_name
                    # bl.write(link_name)
                    # bl.write('\t')
                    self.logger.info(
                        'succeed in linking {} to {}'.format(source, link_name))
                text = '\t' + text1 + '\t' + text2
                bl.write(sample_name)
                bl.write(text)
                bl.write('\n')
            for sn in self.sp_info.keys():
                with open(os.path.join(self.output_dir,"{}_raw_input.txt".format(sn)),"w") as raw:
                    raw.write("RAW\t{}".format("\t".join(self.sp_info[sn])))
            self.logger.info('finish set_output at {}'.format(
                self.__class__.__name__))
        self.end()

    def end(self):
        super(QcModule, self).end()


class TestFunction(unittest.TestCase):
    """
    测试脚本用
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'Qc_' + str(random.randint(1, 10000)),
            "type": "module",
            "name": "tool_lab.ysource.qc",
            "options": {
                "sample_list": "/mnt/ilustre/users/sanger-dev/workspace/20201222/Single_datapre_5046/Datapre/output/sample_list"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
