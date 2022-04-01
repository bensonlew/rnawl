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

class MakelistAgent(Agent):
    def __init__(self, parent):
        super(MakelistAgent, self).__init__(parent)
        options = [
            {"name": "fastq1", "type":"infile","format":"whole_transcriptome.fastq"},
            {"name":"fastq2", "type":"infile","format":"whole_transcriptome.fastq"},
            {"name":"sample_name","type":"string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        if not self.option("fastq1"):
            raise OptionError("请输入样本list")
        if not self.option("fastq2"):
            raise OptionError("请输入样本list")
    
    
    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(MakelistAgent, self).end()

class MakelistTool(Tool):
    def __init__(self, config):
        super(MakelistTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.raw_data = os.path.join(self.software_dir, 'bioinfo/tool_lab/yoogene/raw_data')

    def run(self):
        super(MakelistTool, self).run()
        self.run_makelist()
        self.set_output()
        self.end()

    def run_makelist(self):
        fq1_path = self.option("fastq1").prop["path"]
        fq2_path = self.option("fastq2").prop["path"]
        fq1_name = os.path.basename(fq1_path)
        fq2_name = os.path.basename(fq2_path)
        # self.fq1_sn = fq1_name.rstrip("_1.fq.gz")
        # self.fq2_sn = fq2_name.rstrip("_2.fq.gz")
        # if self.fq1_sn != self.fq2_sn:
        #     self.set_error("fq1和fq2文件不是同一个样本名")
        with open("{}_sample_list".format(self.option("sample_name")),"w") as sl:
            text = "./{}\t{}\t{}\t{}".format(self.option("sample_name"),self.option("sample_name"),fq1_path,fq2_path)
            sl.write(text)

    def set_output(self):
        """
        设置结果目录
        :return:
        """    
        if len(self.output_dir) >0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir, "{}_sample_list".format(self.option("sample_name"))),os.path.join(self.output_dir, "{}_sample_list".format(self.option("sample_name"))))
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
            "id": "makelist_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.makelist",
            "options": dict(
                fastq1="/mnt/ilustre/users/sanger-dev/workspace/20201102/Single_bam2fqgz_7922/Bam2fqgz/output/YGB202000306_1.fq.gz",
                fastq2="/mnt/ilustre/users/sanger-dev/workspace/20201102/Single_bam2fqgz_7922/Bam2fqgz/output/YGB202000306_2.fq.gz",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()