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


class MkbamAgent(Agent):
    def __init__(self, parent):
        super(MkbamAgent, self).__init__(parent)
        options = [
            {"name": "fq1", "type": "infile", "format": "whole_transcriptome.fastq"},
            {"name": "fq2", "type": "infile", "format": "whole_transcriptome.fastq"},
            
            {"name": "sample_name", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        if not self.option("fq1"):
            raise OptionError("因无法找到fq压缩文件，无法生成bam文件")
        if not self.option("fq2"):
            raise OptionError("因无法找到fq压缩文件，无法生成bam文件")
        if not self.option("sample_name"):
            raise OptionError("样本名没有接收到，请检查")

    def set_resource(self):
        self._cpu = 2
        # memory = os.path.getsize(self.option(
        #     'fq1').prop['path'].split(" ")[0])*2
        # self._memory = '{}G'.format(int(memory))
        self._memory = "10G"

    def end(self):
        super(MkbamAgent, self).end()


class MkbamTool(Tool):
    def __init__(self, config):
        super(MkbamTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.samtools_path = "/bioinfo/align/samtools-1.8/samtools"
        self.bwa_path = self.config.SOFTWARE_DIR + "/bioinfo/align/bwa-0.7.17/bwa"
        self.file = {
            'ref': os.path.join(self.software_dir, 'database/Tool_lab/ysource/reference/hg38modMt.fa'),
            'ref_fai': os.path.join(self.software_dir, 'database/Tool_lab/ysource/reference/hg38modMt.fa.fai')
            # 'badpos': os.path.join(self.software_dir, 'database/Tool_lab/ysource/badpos'),
            # 'tree': os.path.join(self.software_dir, 'database/Tool_lab/ysource/tree'),
        }

    def run(self):
        super(MkbamTool, self).run()
        self.run_makesam()
        self.run_sam2bam()
        self.set_output()
        self.end()

    def run_makesam(self):
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        fq1_file = self.option("fq1").prop["path"]
        sample_name = self.option("sample_name")
        fq2_file = self.option("fq2").prop["path"]
        file_path = script_path + "mkbam_{}.sh".format(now_time)
        cmd = "{} mem -t 10 -M -R \"@RG\\tID:{}\\tLB:LB1\\tSM:{}\\tPL:ILLUMINA\" {} {} {} > aln-pe.sam".format(
            self.bwa_path,sample_name,sample_name,self.file["ref"],fq1_file, fq2_file)
        self.logger.info(cmd)
        with open(file_path,"w") as w:
            w.write('#!/bin/bash'+"\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/tool_lab/Yoogene/script_temp/{}".format(
            os.path.basename(file_path))
        self.logger.info("开始生成bam文件")
        self.logger.info(shell)
        command1 =self.add_command("makebamfile", shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("bam生成完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))

        
    
    def run_sam2bam(self):
        sam_file = os.path.join(self.work_dir, "aln-pe.sam")
        cmd = "{} view -bt {} -q 15 -@ 10 {} -o aln-pe.bam".format(self.samtools_path,self.file["ref_fai"],sam_file)
        command1 = self.add_command("samtoolsview", cmd).run()
        self.wait(command1)
        self.logger.info("start samtools views")
        if command1.return_code == 0:
            self.logger.info("bam 生成完成")
        else:
            self.set_error("bam生成出错")
                     
    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir, "aln-pe.bam"),os.path.join(self.output_dir, "aln-pe.bam"))
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
            "id": "mkbam_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.mkbam",
            "options": dict(
                fq1 = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/QC_1.fastq",
                fq2 = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/QC_2.fastq",
                sample_name = "YG202003761",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
