# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.gene_structure.gff3 import Gff3File
from mbio.files.sequence.file_sample import FileSampleFile
from mbio.files.sequence.fasta import FastaFile
from biocluster.config import Config
import pandas as pd
from collections import OrderedDict
import os
import re


class FilecheckItraqAgent(Agent):
    """
    version 1.0
    用于在itraq_and_tmt的workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FilecheckItraqAgent, self).__init__(parent)
        options = [
            {"name": "protein_group", "type": "infile", "format": "itraq_and_tmt.group_table"},
            #分组文件
            {"name": "protein_control", "type": "infile", "format": "itraq_and_tmt.compare_table"},
            #对照组文件
            {"name": "protein_fasta", "type": "infile", 'format': "itraq_and_tmt.common"},
            #蛋白FASTA文件
            {"name": "ratio_exp", "type": "infile", 'format': "itraq_and_tmt.ratio_exp"},
            #蛋白ratio定量表
            {"name": "protein", "type": "infile", 'format': "itraq_and_tmt.common"},
            #蛋白鉴定表
            {"name": "psm", "type": "infile", 'format': "itraq_and_tmt.common"},
            #蛋白PSM表
            {"name": "protein_information", "type": "infile", 'format': "itraq_and_tmt.common"},
            #蛋白信息表
            {"name": "peptide", "type": "infile", 'format': "itraq_and_tmt.common"},
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
            raise OptionError("必须输入蛋白鉴定表", code = "32502801")
        if not self.option('psm'):
            raise OptionError("必须输入蛋白PSM表", code = "32502802")
        if not self.option('protein_information'):
            raise OptionError("必须输入蛋白信息表", code = "32502803")
        if not self.option('peptide'):
            raise OptionError("必须输入肽段表", code = "32502804")
        if not self.option('ratio_exp'):
            raise OptionError("必须输入蛋白Ratio定量表", code = "32502805")
        if not self.option('protein_group'):
            raise OptionError("必须输入分组文件", code = "32502806")
        if not self.option('protein_control'):
            raise OptionError("必须输入对照组文件", code = "32502807")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "5G"


class FilecheckItraqTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(FilecheckItraqTool, self).__init__(config)
        self.ratio_samples = list()
        self.scaled_samples = list()

    def check_exp(self):
        self.logger.info('正在检测ratio定量表')
        self.ratio_samples = self.option("ratio_exp").prop["samples"]
        with open(self.option("ratio_exp").prop["path"], 'r') as exp_r:
            exp = exp_r.read()
            if u'\t\t' in exp:
                self.set_error("蛋白Ratio定量表出错, 里面应该是有空值，另外别忘了改fasta文件", code = "32502808")
        for ss in self.scaled_samples:
            if ss not in self.ratio_samples:
                self.set_error("蛋白Ratio定量表出错, 样本%s在蛋白Ratio定量表文件中未出现", variables = (ss), code = "32502809")
        self.logger.info('检测ratio定量表完毕')

    def check_group(self):
        if self.option('protein_group').is_set:
            self.logger.info("正在检测分组文件")
            self.option("protein_group").get_info()
            gp_sample = self.option("protein_group").prop["sample"]
            for gp in gp_sample:
                if u'-' in gp or type(gp) != str:
                    self.set_error("样品名命名规则：里面不能有‘-’，且不能是数字", code = "32502810")
                if gp not in self.ratio_samples:
                    self.set_error("group表出错, 样本%s在蛋白Ratio定量表中未出现", variables = (gp), code = "32502811")
            gp = self.option("protein_group").prop["group_dict"]
            for g in gp:
                if u'-' in g or type(g) != str:
                    self.set_error("样品组命名规则：里面不能有‘-’，且不能是数字", code = "32502812")
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
                            self.set_error("对照组文件出错，分组%s在分组文件中未出现", variables = (cp), code = "32502813")
        else:
            self.logger.info("未检测到对照组文件， 跳过...")
        self.logger.info("对照组文件检测完毕")

    # 检测protein.xls文件，fasta文件和exp文件，----by 封一统 2018-7-27
    def check_protein_exp_fasta(self):
        protein_list = set()
        exp_list = set()
        fasta_list = set()
        except_list = set()
        with open(self.option('protein').prop['path'], 'r') as protein_r:
            nf = len(protein_r.readline().strip().split('\t'))
            for line in protein_r.readlines():
                line = line.rstrip('\n').split('\t')
                if len(line) > 4:
                    protein_list.add(line[3])
                    if len(line) != nf:
                        self.set_error("protein文件有问题，%s所在的那一行有%s列，而表头有%s列", variables = (line[3], str(len(line)), str(nf)), code = "32502814")
        with open(self.option('protein_fasta').prop['path'], 'r') as fasta_r:
            for line in fasta_r.readlines():
                if line.startswith('>'):
                    fasta_list.add(line.strip().lstrip('>').split(' ')[0])
        with open(self.option('ratio_exp').prop['path'], 'r') as exp_r:
            _ = exp_r.readline()
            for line in exp_r.readlines():
                line = line.strip('\n').split('\t')
                if len(line) > 1:
                    exp_list.add(line[0])
        print(fasta_list==exp_list)
        check1 = exp_list.difference(fasta_list)
        check2 = fasta_list.difference(exp_list)
        if check1:
            self.set_error("有%s个蛋白在exp文件里有，而fasta文件里没有，如%s", variables = (str(len(check1)), list(check1)[0]), code = "32502814")
        if check2:
            self.set_error("有%s个蛋白在fasta文件里有，而exp文件里没有，如%s", variables = (str(len(check2)), list(check2)[0]), code = "32502815")
        for acc in exp_list:
            if not acc in list(protein_list):
                except_list.add(acc)
        if except_list:
            self.set_error("共有%s个蛋白在fasta和exp文件里有，而protein文件里没有，如%s", variables = (str(len(except_list)), list(except_list)[0]), code = "32502816")

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
        # num = df.loc[:, 'Annotated Sequence'].tolist()
        # mylist = [len(x.split(".")[1]) if u'.' in x else len(x) for x in num]
        # mylist.sort()
        # length_dict =OrderedDict()
        # for i in mylist:
        #     length_dict[str(i)] = mylist.count(i)

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
        super(FilecheckItraqTool, self).run()
        self.check_exp()
        self.check_group()
        self.check_control()
        self.check_protein_exp_fasta()
        self.check_psm()
        self.check_peptide()
        self.check_protein_info()
        self.end()
