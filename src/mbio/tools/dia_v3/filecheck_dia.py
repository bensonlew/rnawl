# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
from mbio.files.sequence.fastq import FastqFile
from mbio.files.gene_structure.gff3 import Gff3File
from mbio.files.sequence.file_sample import FileSampleFile
from mbio.files.sequence.fasta import FastaFile
from biocluster.config import Config
import os
import re


class FilecheckDiaAgent(Agent):
    """
    version 1.0
    用于在DIA的workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FilecheckDiaAgent, self).__init__(parent)
        options = [
            {"name": "protein_group", "type": "infile", "format": "labelfree.group_table"},
            #分组文件
            {"name": "protein_control", "type": "infile", "format": "labelfree.compare_table"},
            #对照组文件
            {"name": "protein_fasta", "type": "infile", 'format': "labelfree.common"},
            #蛋白FASTA文件
            {"name": "ratio_exp", "type": "infile", 'format': "labelfree.ratio_exp"},
            #蛋白ratio定量表
            {"name": "protein", "type": "infile", 'format': "labelfree.common"},
            #蛋白鉴定表
            {"name": "psm", "type": "infile", 'format': "labelfree.common"},
            #蛋白PSM表
            {"name": "protein_information", "type": "infile", 'format': "labelfree.common"},
            #蛋白信息表
            {"name": "peptide", "type": "infile", 'format': "labelfree.common"},
            #肽段表
        ]
        self.add_option(options)
        self.step.add_steps("file_check")
        self.on('start', self.start_file_check)
        self.on('end', self.end_file_check)

    def start_file_check(self):
        self.step.file_check.start()
        self.step.update()

    def end_file_check(self):
        self.step.file_check.finish()
        self.step.update()

    def check_option(self):
        if not self.option('protein'):
            raise OptionError("必须输入蛋白鉴定表")
        if not self.option('psm'):
            raise OptionError("必须输入蛋白PSM表")
        if not self.option('protein_information'):
            raise OptionError("必须输入蛋白信息表")
        if not self.option('peptide'):
            raise OptionError("必须输入肽段表")
        if not self.option('ratio_exp'):
            raise OptionError("必须输入蛋白Ratio定量表")
        if not self.option('protein_group'):
            raise OptionError("必须输入分组文件")
        if not self.option('protein_control'):
            raise OptionError("必须输入对照组文件")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "5G"


class FilecheckDiaTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(FilecheckDiaTool, self).__init__(config)
        self.ratio_samples = list()

    def check_exp(self):
        self.logger.info('正在检测ratio定量表')
        self.ratio_samples = self.option("ratio_exp").prop["samples"]
        self.logger.info('检测ratio定量表完毕')

    def check_group(self):
        if self.option('protein_group').is_set:
            self.logger.info("正在检测分组文件")
            self.option("protein_group").get_info()
            gp_sample = self.option("protein_group").prop["sample"]
            for gp in gp_sample:
                if gp not in self.ratio_samples:
                    raise Exception("group表出错, 样本{}在蛋白Ratio定量表中未出现".format(gp))
        else:
            self.logger.info("未检测到分组文件， 跳过...")
        self.logger.info("分组文件检测完毕")

    def check_control(self):
        if self.option('protein_control').is_set:
            self.logger.info("正在检测对照组文件")
            vs_list = self.option("protein_control").prop["cmp_list"]
            con_samples = []
            if self.option('protein_group').is_set:
                group_scheme = self.option('protein_group').prop['group_scheme'][0]
                group_name = self.option('protein_group').get_group_name(group_scheme)
            for vs in vs_list:
                for i in vs:
                    if i not in con_samples:
                        con_samples.append(i)
                for cp in con_samples:
                    if self.option('protein_group').is_set:
                        if cp not in group_name:
                            raise Exception("对照组文件出错，分组{}在分组文件中未出现".format(cp))
        else:
            self.logger.info("未检测到对照组文件， 跳过...")
        self.logger.info("对照组文件检测完毕")

    def run(self):
        super(FilecheckDiaTool, self).run()
        self.check_exp()
        self.check_group()
        self.check_control()
        self.end()
