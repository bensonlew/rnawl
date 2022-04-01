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
import random
import re
import os
import sys
import shutil

class FastqstatAgent(Agent):
    def __init__(self, parent):
        super(FastqstatAgent, self).__init__(parent)
        options = [
            {"name":"raw_input","type":"string"},
            {"name":"clean_input","type":"string"},
            {"name":"sample_name","type":"string"}
        ]
        self.add_option(options)
    
    def check_options(self):
        """
        检查参数
        """
        if not os.path.exists(self.option("raw_input")):
            raise OptionError("找不到raw_input.txt文件")
        if not os.path.exists(self.option("clean_input")):
            raise OptionError("找不到clean_input文件")
        return True
    
    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'
    
    def end(self):
        super(FastqstatAgent, self).end()

class FastqstatTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(FastqstatTool, self).__init__(config)
        self.software_dir = self.config.SOFTWARE_DIR
        self.Fastqstat_path = self.config.SOFTWARE_DIR + "/bioinfo/seq/FastqStat.jar"
        self.java_path = self.config.SOFTWARE_DIR + "/program/sun_jdk1.8.0/bin/"
        self.raw_input = self.option("raw_input")
        self.clean_input = self.option("clean_input")

    def fastqstat(self,infile,outfile):
        if not os.path.exists("./QC_1.fastq"):
            os.link("/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/QC_1.fastq","./QC_1.fastq")
        if not os.path.exists("./QC_2.fastq"):
            os.link("/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/QC_1.fastq","./QC_2.fastq")
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "fastqstat_{}.sh".format(now_time)
        cmd = "{}java -jar {} -i {} > {}_statistic.xls".format(self.java_path,self.Fastqstat_path,infile,outfile)
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
        self.logger.info("开始生成fastq统计文件")
        self.logger.info(shell)
        command1 =self.add_command("fastqstat{}".format(outfile), shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行fastqstat完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))
        # command = self.add_command("fastqsstat", cmd, ignore_error = True)
        # command.run()
        # self.wait(command)
        # if command.return_code == 0:
        #     self.logger.info("运行Fastqstat.jar完成")
        #     os.system("awk -F '\t' BEGIN'{OFS=\"\t\"}''{print \$1,\$2,\$3,\$11,\$12,\$13,\$14}' {}_statistic.xls > ./{}_statistic_reform.xls".format(outfile, outfile))
        #     # self.set_output()
        # elif command.return_code in [1, -9]:
        #     self.add_state("memory_limit", "memory is low!")
        # else:
        #     self.set_error("运行Fastqstat.jar运行出错！")
        #     return False

    
    def get_reform(self, inputt):
        cmd = "awk -F '\\t' BEGIN'{OFS=\"\\t\"}''{print $1,$2,$3,$11,$12,$13,$14}'"
        cmd = cmd + " {}_statistic.xls > ./{}_statistic_reform.xls".format(inputt, inputt)
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "reform_{}.sh".format(now_time)
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
        self.logger.info("开始生成fastq统计文件")
        self.logger.info(shell)
        command1 =self.add_command("reform{}".format(inputt.lower()), shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行fastqstat完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))

        # command = self.add_command("awk{}".format(input), cmd, ignore_error = True)
        # command.run()
        # self.wait(command)
        # if command.return_code == 0:
        #     self.logger.info("{}_statistic_reform.xls创建成功".format(input))
        #     # self.set_output()
        # else:
        #     self.set_error("{}_statistic_reform.xls创建失败".format(input))
        #     return False

    def set_output(self):
        self.logger.info("set output")
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir, "raw_statistic_reform.xls"),
                    os.path.join(self.output_dir, "{}_raw_statistic_reform.xls".format(self.option('sample_name'))))
            os.link(os.path.join(self.work_dir, "clean_statistic_reform.xls"),
                    os.path.join(self.output_dir, "{}_clean_statistic_reform.xls".format(self.option('sample_name'))))
            # os.link(os.path.join(self.work_dir, "raw_statistic.xls"),
            #         os.path.join(self.output_dir, "raw_statistic.xls"))
            # os.link(os.path.join(self.work_dir, "clean_statistic.xls"),
            #         os.path.join(self.output_dir, "clean_statistic.xls"))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))
        
    def run(self):
        super(FastqstatTool, self).run()
        self.fastqstat(self.option("raw_input"),"raw")
        self.fastqstat(self.option("clean_input"),"clean")
        self.get_reform("raw")
        self.get_reform("clean")
        self.set_output()
        self.end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """
        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            data = {
                "id": "fastqstat_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "tool_lab.ysource.fastqstat",
                "instant": False,
                "options": {
                    "raw_input":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/raw_input.txt",
                    "clean_input":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/clean_input.txt"
                }
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()



    

