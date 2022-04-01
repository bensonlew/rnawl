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

class Bam2fqgzAgent(Agent):
    def __init__(self, parent):
        super(Bam2fqgzAgent, self).__init__(parent)
        options= [
            {"name":"bam_file", "type":"infile","format":"align.bwa.bam"},
            {"name":"sample_name","type":"string"}
        ]
        self.add_option(options)

    def check_option(self):
        '''
        参数检查
        '''
        if not self.option("bam_file"):
            raise OptionError("请选择需要转换的bam文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'
    
    def end(self):
        super(Bam2fqgzAgent, self).end()

class Bam2fqgzTool(Tool):
    def __init__(self, config):
        super(Bam2fqgzTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.samtools = os.path.join(self.software_dir, 'bioinfo/align/samtools-1.8')
        self.bam = self.option("bam_file").prop['path']
        self.bam_name = self.option("sample_name")
        # self.bam_name = os.path.basename(self.option("bam_file"))
        
    def run(self):
        super(Bam2fqgzTool, self).run()
        self.run_bam2fqgz()
        self.set_output()
        self.end()
        
    def run_bam2fqgz(self):
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "bam2fq_{}.sh".format(now_time)
        cmd = "{} sort -n {} | {} fastq -1 {}_1.fq.gz -2 {}_2.fq.gz -s {}.fq -".format(self.samtools,self.bam,self.samtools,self.bam_name,self.bam_name,self.bam_name)
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash'+"\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/tool_lab/Yoogene/script_temp/{}".format(
            os.path.basename(file_path))
        self.logger.info("开始生成fastq压缩文件")
        self.logger.info(shell)
        command1 =self.add_command("bamtofqgz", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行samtools完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir, "{}_1.fq.gz".format(self.bam_name)),
                    os.path.join(self.output_dir, "{}_1.fq.gz".format(self.bam_name)))
            os.link(os.path.join(self.work_dir, "{}_2.fq.gz".format(self.bam_name)),
                    os.path.join(self.output_dir, "{}_2.fq.gz".format(self.bam_name)))
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
            "id": "bam2fqgz_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.bam2fqgz",
            "options": {
                "bam_file":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YGB202000306.bam",
                "sample_name":"YGB202000306"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()