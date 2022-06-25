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

class MappingAgent(Agent):
    def __init__(self, parent):
        super(MappingAgent, self).__init__(parent)
        options= [
            {"name":"bam_file","type":"infile","format":"align.bwa.bam"},
            {"name":"sample_name","type":"string"}
        ]
        self.add_option(options)

    def check_option(self):
        '''
         参数检查
        '''
        if not self.option("bam_file"):
            raise OptionError("没有找到bam文件，请检查上一步")

    def set_resource(self):
        self._cpu = 2
        self._memory = '30G'
    
    def end(self):
        super(MappingAgent, self).end()

class MappingTool(Tool):
    def __init__(self, config):
        super(MappingTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.samtools = 'miniconda2/bin/samtools'
        self.bam = self.option("bam_file").prop["path"]

    def run(self):
        super(MappingTool, self).run()
        self.run_samtoolsort()
        self.run_samtoolrmdup()
        self.run_samtoolindex()
        self.set_output()
        self.end()

    def run_samtoolsort(self):
        self.sort_bam = "aln-pe.sort.bam"
        samtools = os.path.join(self.software_dir, 'miniconda2/bin/samtools')
        cmd = "{} sort -@ 10 {} >{}".format(samtools,self.bam,self.sort_bam)
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.software_dir + '/bioinfo/tool_lab/Yoogene/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "samtoolssort_{}.sh".format(now_time)
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
        command1 =self.add_command("samtoolssort{}".format(self.option('sample_name').lower()), shell).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行samtools完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")
        os.system('rm {}'.format(file_path))

    
    def run_samtoolrmdup(self):
        self.rmdup_bam = "align_final_sort.bam"
        cmd = "{} rmdup {} {}".format(self.samtools, self.sort_bam,self.rmdup_bam)
        command1 = self.add_command("samtoolrmdup{}".format(self.option('sample_name').lower()),cmd).run()
        self.wait(command1)
        self.logger.info("start samtools rmdup")
        if command1.return_code == 0:
            self.logger.info("bam去重复完成")
        else:
            self.set_error("samtools rmdup 出错")
    
    def run_samtoolindex(self):
        cmd = "{} index {}".format(self.samtools, self.rmdup_bam)
        command1 = self.add_command("samtoolindex{}".format(self.option('sample_name').lower()), cmd).run()
        self.wait(command1)
        self.logger.info("start samtools index")
        if command1.return_code == 0:
            self.logger.info("bam建立索引完成")
        else:
            self.set_error("samtools index 出错")
        cmd1 = "{} index {}".format(self.samtools, self.sort_bam)
        command2 = self.add_command("samtoolindex2_{}".format(self.option('sample_name').lower()), cmd1).run()
        self.wait(command2)
        self.logger.info("start samtools index")
        if command2.return_code == 0:
            self.logger.info("bam建立索引完成")
        else:
            self.set_error("samtools index 出错")
    
    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(os.path.join(self.work_dir,"aln-pe.sort.bam"),os.path.join(self.output_dir,"aln-pe.sort.bam"))
            os.link(os.path.join(self.work_dir,"aln-pe.sort.bam.bai"),os.path.join(self.output_dir,"aln-pe.sort.bam.bai"))
            os.link(os.path.join(self.work_dir,"align_final_sort.bam.bai"),os.path.join(self.output_dir,"align_final_sort.bam.bai"))
            os.link(os.path.join(self.work_dir,"align_final_sort.bam"),os.path.join(self.output_dir,"align_final_sort.bam"))
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
            "id": "mapping_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.mapping",
            "options": {
                "bam_file":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/QC_and_mapping_temp/aln-pe.bam",
                "sample_name":"YG202003761"
                }
            }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()