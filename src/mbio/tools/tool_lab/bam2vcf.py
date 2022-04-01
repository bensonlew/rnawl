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


class Bam2vcfAgent(Agent):
    def __init__(self, parent):
        super(Bam2vcfAgent, self).__init__(parent)
        options = [
            # {"name":"sample1","type":"infile","format":"align.bwa.bam"},#用作共祖时间计算的样本1，要提供单倍型
            # {"name":"sample2","type":"infile","format":"align.bwa.bam"},#用作共祖时间计算的样本2
            # {"name":"outgroup1","type":"infile","format":"align.bwa.bam"},#用作外群的bam文件
            # {"name":"outgroup2","type":"infile","format":"align.bwa.bam"},
            # {"name":"outgroup3","type":"infile","format":"align.bwa.bam"},
            # {"name":"outgroup4","type":"infile","format":"align.bwa.bam"},
            # {"name":"outgroup5","type":"infile","format":"align.bwa.bam"},
            {"name": "og_bam_dir", "type": "infile",
                "format": "align.bwa.bam_dir"},  # 存放外群的bam文件的文件夹
            {"name": "sample1", "type": "infile",
                "format": "align.bwa.bam"},  # 用作共祖时间计算的样本1，要提供单倍型
            {"name": "c_bam_dir", "type": "infile",
                "format": "align.bwa.bam_dir"},  # 存放用于比较的样本的bam文件的文件夹
        ]
        self.add_option(options)

    def check_options(self):
        '''
        参数检查
        '''
        if not self.option("sample1"):
            raise OptionError("请选择需要计算共祖时间的样品1")
        if not self.option("og_bam_dir"):
            raise OptionError("请选择用于计算与sample1共祖时间的样品所在的文件夹")
        if not self.option("c_bam_dir"):
            raise OptionError("请选择外群样品所在的文件夹")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(Bam2vcfAgent, self).end()


class Bam2vcfTool(Tool):
    def __init__(self, config):
        super(Bam2vcfTool, self).__init__(config)
        self._version = "v1.0"
        self.software_dir = self.config.SOFTWARE_DIR
        self.samtools = os.path.join(
            self.software_dir, 'program/Python/bin/samtools')
        self.bcftools = os.path.join(
            self.software_dir, 'bioinfo/align/bcftools-1.6/bcftools')
        self.file = {
            'ref_fa': os.path.join(self.software_dir, 'database/Tool_lab/ysource/ref.fa'),
            'bed': os.path.join(self.software_dir, 'database/Tool_lab/ysource/hg38_8_43M.bed'),

        }

    def run(self):
        super(Bam2vcfTool, self).run()
        self.run_mvbamfile()
        self.run_bam2vcf()
        self.set_output()
        self.end()

    def run_mvbamfile(self):
        self.logger.info("开始创建bam.list文件")
        # self.bamdir_path = self.option("bam_dir").prop["path"]
        '''
        创建一个bam.list
        '''
        s1 = "{}".format(self.option("sample1").prop["path"])
        self.logger.info("sample: " + s1)
        bamlist = s1
        bamlist1 = "{}".format(os.path.basename(
            self.option("sample1").prop["path"]))
        # s2 = "{}".format(self.option("sample2").prop["path"])
        # t1 = "{}".format(os.path.basename(self.option("sample1").prop["path"]))
        # t2 = "{}".format(os.path.basename(self.option("sample1").prop["path"]))
        c_bam_dir = "{}".format(self.option("c_bam_dir").prop["path"])
        for root, _dir, files in os.walk(c_bam_dir):
            for name in files:
                self.logger.info("增加比较样本：{}".format(name))
                bamlist += "\n{}".format(os.path.join(root, name))
                bamlist1 += "\n{}".format(name)
        self.logger.info("add sample used to count:" + bamlist)
        outgroup_dir = "{}".format(self.option("og_bam_dir").prop["path"])
        for root, _dir, files in os.walk(outgroup_dir):
            for name in files:
                self.logger.info("增加外群：{}".format(name))
                bamlist += "\n{}".format(os.path.join(root, name))
                bamlist1 += "\n{}".format(name)
        self.logger.info("bamlist加上外群后：" + bamlist)
        # self.out_group = ""
        # self.out_group = '{}'.format(outgroup1)

        # self.logger.info("{}".format(self.out_group))
        self.bam_list = open('bam.list', 'w')
        self.bam_list1 = open('bam1.list', 'w')
        content = bamlist

        content1 = bamlist1
        self.bam_list.write(content)
        self.bam_list1.write(content1)
        self.logger.info("写入bam.list文件")
        self.logger.info(content)
        self.logger.info(content1)
        self.bam_list.close()
        self.bam_list1.close()
        if not os.path.exists(os.path.join(self.work_dir, 'bam.list')):
            self.set_error("bam.list没有生成")
        # self.bamlist = self.out_group.rstrip().split(',')
        # self.bamlist.append(s1)
        # self.bamlist.append(s2)

        # with open(os.path.join(self.work_dir, 'bam.list')) as bl:
        #     for line in bl:
        #         if line.lstrip('{}'.format(self.bamdir_path)).rstrip('\n') not in self.bamlist:
        #             self.set_error("错误，bamlist写入失败")
        #         else:
        #             self.logger.info("写入成功")

    def run_bam2vcf(self):
        now_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f") + \
            "_" + str(random.randint(1, 10000))
        script_path = self.config.SOFTWARE_DIR + '/bioinfo/ancestorage/script_temp/'
        if not os.path.exists(script_path):
            os.mkdir(script_path)
        file_path = script_path + "bam2vcf_{}.sh".format(now_time)
        cmd = "{} mpileup -t AD,DP -b {} -uf {} -l {}|{} call -mv > {}".format(
            self.samtools, os.path.join(self.work_dir, 'bam.list'), self.file['ref_fa'], self.file['bed'], self.bcftools, 'samples-filter-8.4.vcf')
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
        command1 = self.add_command("bam2vcf", shell).run()
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
            os.link(os.path.join(self.work_dir, "samples-filter-8.4.vcf"),
                    os.path.join(self.output_dir, "samples-filter-8.4.vcf"))
            os.link(os.path.join(self.work_dir, "bam.list"),
                    os.path.join(self.output_dir, "bam.list"))
            os.link(os.path.join(self.work_dir, "bam1.list"),
                    os.path.join(self.output_dir, "bam1.list"))
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
            "id": "bam2vcf_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.bam2vcf",
            "options": dict(

                # bam_dir="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/de_tools/bam_dir",
                sample1="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/de_tools/bam_dir/YG201902325.hg38.chrY.dedup.sort.bam",
                og_bam_dir="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/de_tools/outgroup",
                c_bam_dir="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/de_tools/c_bam",
                
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
