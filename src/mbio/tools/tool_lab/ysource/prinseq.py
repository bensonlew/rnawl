# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from collections import namedtuple, defaultdict
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import random
import datetime
import subprocess
import re
import os
import sys
import shutil

class PrinseqAgent(Agent):
    def __init__(self, parent):
        super(PrinseqAgent, self).__init__(parent)
        options = [
            # {"name":"f1","type":"infile","format":"sequence.fastq"},
            # {"name":"f2","type":"infile","format":"sequence.fastq"},
            {"name":"f1","type":"string"},
            {"name":"f2","type":"string"},
            {"name":"sample_name","type":"string"},
            {"name":"quality_score","type":"int","default":20},
            {"name":"min_length","type":"int","default": 25},
            {"name":"trim_ns","type":"int","default": 1},
            {"name":"trim_left","type":"int","default": 0},
            {"name":"ns_max_n","type":"int","default":5},
            {"name":"trim_qual_window","type":"int","default":5},
            {"name":"trim_qual_step","type":"int","default":1},
            {"name":"line_width","type":"int","default":1}
        ]
        self.add_option(options)
    
    def check_options(self):
        '''
        参数检查
        '''
        if not os.path.exists(self.option("f1")):
            raise OptionError("请输入原始文件")
        if not os.path.exists(self.option("f2")):
            raise OptionError("请输入原始文件")

    def set_resource(self):
        """
        设置所需资源
        """    
        self._cpu = 2
        self._memory = '30G'
    
    def end(self):
        super(PrinseqAgent, self).end()

class PrinseqTool(Tool):
    def __init__(self, config):
        super(PrinseqTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.prinseq_path = self.software_dir + "/bioinfo/medical/prinseq-lite-0.20.4"
        self.Prinseq = os.path.join(self.prinseq_path,'prinseq-lite.pl')
        


    
    def run(self):
        super(PrinseqTool, self).run()
        self.run_Prinseq()
        self.set_output()
        self.end()

    def run_Prinseq(self):
        cmd = "{} ".format(self.Prinseq)
        self.path1 = self.option('f1')
        self.path2 = self.option('f2')
        cmd += "-fastq {} -fastq2 {} ".format(self.path1, self.path2)
        trim_ns = self.option("trim_ns")
        trim_left =self.option("trim_left")
        ns_max_n = self.option("ns_max_n")
        cmd += "-trim_ns_left {} -trim_ns_right {} -trim_left {} -ns_max_n {} ".format(trim_ns,trim_ns,trim_left,ns_max_n)
        quality_score = self.option('quality_score')
        min_length = self.option('min_length')
        cmd += "-trim_qual_right {} -trim_qual_left {} -min_len {} ".format(quality_score, quality_score, min_length)
        trim_qual_window = self.option("trim_qual_window")
        trim_qual_step = self.option("trim_qual_step")
        line_width = self.option("line_width")
        cmd += "-trim_qual_window {} -trim_qual_step {} -out_good {}_QC _out_bad null -line_width {} -no_qual_header >> ./Miseq_QC5.log".format(trim_qual_window,trim_qual_step,self.option('sample_name'),line_width)
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "prinseq_{}.sh".format(now_time)
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
        command1 =self.add_command("prinseq", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行prinseq完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))
        # command1 = self.add_command("prinseq", cmd).run()
        # self.wait(command1)
        # self.logger.info("start prinseq")
        # if command1.return_code == 0:
        #     self.logger.info("Prinseq运行完成")
        # else:
        #     self.set_error("Prinseq运行出错!")



    def set_output(self):
        """
        设置结果目录
        :return:
        """
        output_name = "{}_QC".format(self.option('sample_name'))
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        self.logger.info(os.path.join(self.work_dir,"{}_1.fastq".format(output_name)))
        try:
             os.link(os.path.join(self.work_dir,"{}_1.fastq".format(output_name)),os.path.join(self.output_dir,"{}_1.fastq".format(output_name)))
             os.link(os.path.join(self.work_dir,"{}_2.fastq".format(output_name)),os.path.join(self.output_dir,"{}_2.fastq".format(output_name)))
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
            "id": "prinseq_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.prinseq",
            "options": dict(
                f1="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/R1_cutadapt",
                f2="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/R2_cutadapt",
                sample_name="YG202003761"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
