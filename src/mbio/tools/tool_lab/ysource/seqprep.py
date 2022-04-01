# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
# from biocluster.config import config
from collections import namedtuple, defaultdict
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import datetime
import subprocess
import re
import os
import sys
import shutil

class SeqprepAgent(Agent):
    def __init__(self, parent):
        super(SeqprepAgent, self).__init__(parent)
        options = [
            # {"name":"R1","type":"infile","format":"sequence.fastq"},
            # {"name":"R2","type":"infile","format":"sequence.fastq"},
            {"name":"R1","type":"string"},
            {"name":"R2","type":"string"},
            {"name":"sample_name","type":"string"},
            {"name":"quality_score","type":"int","default":20},
            {"name":"min_length","type":"int","default": 25},
            {"name":"f_adapter","type":"string","default":"GATCGGAAGAGCACACGTCT"},
            {"name":"r_adapter","type":"string","default":"AATGATACGGCGACCACCGA"}
        ]
        self.add_option(options)
    
    def check_options(self):
        '''
        参数检查
        '''
        # if not self.option("R1").is_set:
        #     raise OptionError("请输入原始文件")
        # if not self.option("R2").is_set:
        #     raise OptionError("请输入原始文件")
        if not self.option("R1"):
            raise OptionError("请输入原始文件")
        if not self.option("R2"):
            raise OptionError("请输入原始文件")


    def set_resource(self):
        """
        设置所需资源
        """    
        self._cpu = 2
        # memory = os.path.getsize(self.option('R1').prop["path"].split(" ")[0])*2
        self._memory = '30G'
    
    def end(self):
        super(SeqprepAgent, self).end()

class SeqprepTool(Tool):
    def __init__(self, config):
        super(SeqprepTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.seqprep_path = '/bioinfo/seq'
        self.seqprep = os.path.join(self.seqprep_path, 'SeqPrep')
        


    
    def run(self):
        super(SeqprepTool, self).run()
        self.run_seqprep()
        self.run_gzip(1)
        self.run_gzip(2)
        self.set_output()
        self.end()

    def run_seqprep(self):
        cmd = "{} ".format(self.seqprep)
        # self.path1 = self.option('R1').prop["path"]
        # self.path2 = self.option('R2').prop["path"]
        self.path1 = self.option('R1')
        self.path2 = self.option('R2')
        cmd += "-f {} -r {} ".format(self.path1, self.path2)
        cmd += "-1 R1_cutadapt_{}.gz -2 R2_cutadapt_{}.gz ".format(self.option('sample_name'),self.option('sample_name'))
        quality_score = self.option('quality_score')
        min_length = self.option('min_length')
        cmd += "-q {} -L {} ".format(quality_score, min_length)
        f_adapter = self.option('f_adapter')
        r_adapter = self.option('r_adapter')
        cmd += "-A '{}' -B '{}'".format(f_adapter,r_adapter)
        command1 = self.add_command("seqprep", cmd).run()
        self.logger.info("start seqprep")
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("seqprep运行完成")
        else:
            self.set_error("seqprep运行出错!")

    def run_gzip(self,index):
        cmd = "gzip -df R{}_cutadapt_{}.gz".format(index,self.option('sample_name'))
        command1 = self.add_command("gzip{}".format(index), cmd).run()
        self.logger.info("start gzip R{}_cutadapt_{}.gz".format(index,self.option('sample_name')))
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("gzip运行完成")
        else:
            self.set_error("gzip运行出错!")
        


    def set_output(self):
        """
        设置结果目录
        :return:
        """
        output_name = "_cutadapt_{}".format(self.option('sample_name'))
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
             os.link(os.path.join(self.work_dir,"R1{}".format(output_name)),os.path.join(self.output_dir,"R1{}".format(output_name)))
             os.link(os.path.join(self.work_dir,"R2{}".format(output_name)),os.path.join(self.output_dir,"R2{}".format(output_name)))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))
            
class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "seqprep_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.seqprep",
            "options": dict(
                R1="/mnt/ilustre/users/sanger-dev/workspace/20201107/Single_datapre_2967/Datapre/Bam2fqgz/output/YGB202000308_1.fq.gz",
                R2="/mnt/ilustre/users/sanger-dev/workspace/20201107/Single_datapre_2967/Datapre/Bam2fqgz/output/YGB202000308_2.fq.gz",
                sample_name="YGB202000308"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()