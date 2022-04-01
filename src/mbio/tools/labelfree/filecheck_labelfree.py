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


class FilecheckLabelfreeAgent(Agent):
    """
    version 1.0
    用于在labelfree的workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FilecheckLabelfreeAgent, self).__init__(parent)
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


class FilecheckLabelfreeTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(FilecheckLabelfreeTool, self).__init__(config)
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

    def check_psm(self):
        self.logger.info('正在检测PSM表')
        try:
            df = pd.read_table(self.option('psm').prop['path'], sep="\t", header=0)
        except:
            self.set_error('上传的PSM表无法正常打开，请检查文件格式。')
        columns = df.columns.tolist()
        checking_list = ['m/z [Da]', 'DeltaM [ppm]']
        excluded_list = [i for i in checking_list if i not in columns]
        if excluded_list:
            self.set_error("请检查上传的PSM表，已知缺少{}列。".format(','.join(excluded_list)))
        self.logger.info('PSM表检测完毕')

    def check_peptide(self):
        self.logger.info('正在检测肽段表')
        try:
            df = pd.read_table(self.option('peptide').prop['path'], sep ="\t", header = 0)
        except:
            self.set_error('上传的肽段表无法正常打开，请检查文件格式。')
        columns = df.columns.tolist()
        checking_list = ['Annotated Sequence']
        excluded_list = [i for i in checking_list if i not in columns]
        if excluded_list:
            self.set_error("请检查上传的肽段表，已知缺少{}列。".format(','.join(excluded_list)))
        self.logger.info('肽段表检测完毕')

    def check_protein_info(self):
        self.logger.info('正在检测蛋白信息表')
        try:
            df = pd.read_table(self.option('protein_information').prop['path'], sep="\t", header=0)
        except:
            self.set_error('上传的蛋白信息表无法正常打开，请检查文件格式。')
        columns = df.columns.tolist()
        checking_list = [['Identified Spectrum', 'Peptide-Spectrum Matches'],
                         ['Peptide number', 'Peptide sequences'], ['Protein number', 'Proteins'],
                         ['Protein group number', 'Protein groups']]
        excluded_list = [i[0] for i in checking_list if i[0] not in columns and i[1] not in columns]
        if 'Total Spectrum' not in columns:
            excluded_list.append('Total Spectrum')
        if excluded_list:
            self.set_error("请检查上传的蛋白信息表，已知缺少{}列。".format(','.join(excluded_list)))
        self.logger.info('蛋白信息表检测完毕')
        # # 感觉有一次修改之后，列名就对不上了，没办法就只能暴力解决了
        # old_columns = ['Total Spectrum', 'Identified Spectrum', 'Peptide number',
        #                'Protein number', 'Protein group number',
        #                'Peptide-Spectrum Matches', 'Peptide sequences', 'Proteins', 'Protein groups']
        # new_columns = ['total_spectrum', 'identified_spectrum', 'peptide_num', 'protein_num', 'protein_group_num',
        #                'identified_spectrum', 'peptide_num', 'protein_num', 'protein_group_num']
        # columns = zip(old_columns, new_columns)
        # df_col = df.columns.tolist()
        # for o, n in columns:
        #     if o in df_col:
        #         df[n] = df[o]
        # # 咱也不知道前端那边咋弄的啊，还是多加几列吧
        # try:
        #     df['Peptide sequence'] = df['peptide_num']
        #     df['Peptide sequences'] = df['peptide_num']
        #     df['Protein groups'] = df['protein_group_num']
        #     df['Proteins'] = df['protein_num']
        # except:
        #     pass

    def run(self):
        super(FilecheckLabelfreeTool, self).run()
        self.check_exp()
        self.check_group()
        self.check_control()
        self.check_psm()
        self.check_peptide()
        self.check_protein_info()
        self.end()
