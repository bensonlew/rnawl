# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import datetime
import random
import re
import os
import sys
import shutil

class MkchrybamAgent(Agent):
    def __init__(self, parent):
        super(MkchrybamAgent, self).__init__(parent)
        options = [
            {"name":"bam_file","type":"infile",'format':"align.bwa.bam"},
            {"name":"sample_name","type":"string"}
        ]
        self.add_option(options)
    
    def check_option(self):
        '''
        参数检查
        '''
        if not self.option("bam_file"):
            raise OptionError("没有找到bam文件，请检查上一步")
        if not self.option("sample_name"):
            raise OptionError("没有输入样本名")
    
    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"
    
    def end(self):
        super(MkchrybamAgent, self).end()

class MkchrybamTool(Tool):
    def __init__(self, config):
        super(MkchrybamTool,self).__init__(config)
        self._version = "v1.0"
        self.software_dir=self.config.SOFTWARE_DIR
        self.samtools = 'program/Python/bin/samtools'
        self.bam = self.option("bam_file").prop["path"]
        self.sample_name = self.option("sample_name")

    def run(self):
        super(MkchrybamTool, self).run()
        self.run_samview()
        self.run_samsort()
        self.set_output()
        self.end()
        

    def run_samview(self):
        cmd = "{} view -b -@ 10 {} chrY -o chrY.dedup.bam".format(self.samtools,self.bam)
        command1 = self.add_command("samtoolview",cmd).run()
        self.wait(command1)
        self.logger.info("start samtools view for chrY.bam")
        if command1.return_code == 0:
            self.logger.info("samtools view 完成")
        else:
            self.set_error("samtools view 出错")


    def run_samsort(self):
        bam_chry = "chrY.dedup.bam"
        cmd = "{} sort -@ 10 {} -o {}.hg38.chrY.dedup.sort.bam".format(self.samtools,bam_chry,self.sample_name)
        command1 = self.add_command("samtoolsort",cmd).run()
        self.wait(command1)
        self.logger.info("start samtools sort for chrY.bam")
        if command1.return_code == 0:
            self.logger.info("samtools sort 完成")
        else:
            self.set_error("samtools sort 出错")

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir,"{}.hg38.chrY.dedup.sort.bam".format(self.sample_name)),os.path.join(self.output_dir,"{}.hg38.chrY.dedup.sort.bam".format(self.sample_name)))
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
            "id": "mkchrybam_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.mkchrybam",
            "options": dict(
                bam_file = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/2_bwa/align_final_sort.bam",
                sample_name = "YG202003761"
                
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()