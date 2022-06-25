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


class BamtovcfAgent(Agent):
    def __init__(self, parent):
        super(BamtovcfAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile",
                "format": "align.bwa.bam"},
            {"name":"bed","type":"string"},
            {"name":"output_name","type":"string"},
        ]
        self.add_option(options)

    def check_options(self):
        '''
        参数检查
        '''
        if not self.option("bam_file"):
            raise OptionError("缺少输入参数bam文件")
        if not self.option("output_name"):
            raise OptionError("缺少输入参数输出文件名")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(BamtovcfAgent, self).end()


class BamtovcfTool(Tool):
    def __init__(self, config):
        super(BamtovcfTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.samtools = os.path.join(
            self.software_dir, 'miniconda2/bin/samtools')
        self.bcftools = os.path.join(
            self.software_dir, 'bioinfo/align/bcftools-1.6/bcftools')
        self.vcf_name = self.option("output_name")
        self.file = {
            'ref_fa': os.path.join(self.software_dir, 'database/Tool_lab/ysource/ref.fa'),
            'bed': os.path.join(self.software_dir,self.option("bed")),
        }

    def run(self):
        super(BamtovcfTool, self).run()
        # self.run_mvbamfile()
        self.run_bamtovcf()
        self.set_output()
        self.end()

    # def run_mvbamfile(self):
    #     self.logger.info("开始创建bam.list文件")
    #     # self.bamdir_path = self.option("bam_dir").prop["path"]
    #     '''
    #     创建一个bam.list
    #     '''
    #     s1 = "{}".format(self.option("sample1").prop["path"])
    #     self.logger.info("sample: " + s1)
    #     bamlist = s1
    #     bamlist1 = "{}".format(os.path.basename(
    #         self.option("sample1").prop["path"]))
    #     # s2 = "{}".format(self.option("sample2").prop["path"])
    #     # t1 = "{}".format(os.path.basename(self.option("sample1").prop["path"]))
    #     # t2 = "{}".format(os.path.basename(self.option("sample1").prop["path"]))
    #     c_bam_dir = "{}".format(self.option("c_bam_dir").prop["path"])
    #     for root, _dir, files in os.walk(c_bam_dir):
    #         for name in files:
    #             self.logger.info("增加比较样本：{}".format(name))
    #             bamlist += "\n{}".format(os.path.join(root, name))
    #             bamlist1 += "\n{}".format(name)
    #     self.logger.info("add sample used to count:" + bamlist)
    #     outgroup_dir = "{}".format(self.option("og_bam_dir").prop["path"])
    #     for root, _dir, files in os.walk(outgroup_dir):
    #         for name in files:
    #             self.logger.info("增加外群：{}".format(name))
    #             bamlist += "\n{}".format(os.path.join(root, name))
    #             bamlist1 += "\n{}".format(name)
    #     self.logger.info("bamlist加上外群后：" + bamlist)
    #     # self.out_group = ""
    #     # self.out_group = '{}'.format(outgroup1)

    #     # self.logger.info("{}".format(self.out_group))
    #     self.bam_list = open('bam.list', 'w')
    #     self.bam_list1 = open('bam1.list', 'w')
    #     content = bamlist

    #     content1 = bamlist1
    #     self.bam_list.write(content)
    #     self.bam_list1.write(content1)
    #     self.logger.info("写入bam.list文件")
    #     self.logger.info(content)
    #     self.logger.info(content1)
    #     self.bam_list.close()
    #     self.bam_list1.close()
    #     if not os.path.exists(os.path.join(self.work_dir, 'bam.list')):
    #         self.set_error("bam.list没有生成")
    #     # self.bamlist = self.out_group.rstrip().split(',')
    #     # self.bamlist.append(s1)
    #     # self.bamlist.append(s2)

    #     # with open(os.path.join(self.work_dir, 'bam.list')) as bl:
    #     #     for line in bl:
    #     #         if line.lstrip('{}'.format(self.bamdir_path)).rstrip('\n') not in self.bamlist:
    #     #             self.set_error("错误，bamlist写入失败")
    #     #         else:
    #     #             self.logger.info("写入成功")

    def run_bamtovcf(self):
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/ancestorage/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "bamtovcf_{}.sh".format(now_time)
        cmd = "{} mpileup -d 30000  -uf {} {} -l {}|{} call -mA --skip-variants indels  > {}.vcf".format(
            self.samtools, self.file['ref_fa'],self.option("bam_file").prop["path"], self.file['bed'], self.bcftools, self.option("output_name"))
        self.logger.info(cmd)
        with open(file_path, 'w') as w:
            w.write('#!/bin/bash'+"\n")
            w.write(cmd)
        code = os.system('chmod +x {}'.format(file_path))
        if code == 0:
            self.logger.info("修改{}为可执行文件成功！".format(file_path))
        else:
            self.set_error("修改{}为可执行文件失败！".format(file_path))
        shell = "/bioinfo/ancestorage/script_temp/{}".format(
            os.path.basename(file_path))
        self.logger.info("开始生成vcf文件")
        self.logger.info(shell)
        command1 = self.add_command("bamtovcf", shell).run()
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
            os.link(os.path.join(self.work_dir, "{}.vcf".format(self.vcf_name)),
                    os.path.join(self.output_dir, "{}.vcf".format(self.vcf_name)))

        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunctionf(unittest.TestCase):
    '''
    测试脚本
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "bamtovcf_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ysource.bamtovcf",
            "options": dict(
                bam_file= "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/yfull/test_data/YG202003761/2_bwa/align_final_sort.bam",
                output_name="YG202003761",
                bed = "/mnt/ilustre/users/sanger-dev/app/database/Tool_lab/ysource/hg38_8_43M.bed"

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
