# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.core.exceptions import FileError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.gene_structure.gff3 import Gff3File
from mbio.files.sequence.file_sample import FileSampleFile
from mbio.files.sequence.fasta import FastaFile
from biocluster.config import Config
import os
import re


class FilecheckDenovoAgent(Agent):
    """
    version 1.0
    用于在denovorna的workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(FilecheckDenovoAgent, self).__init__(parent)
        options = [
            {"name": "sample_num", "type": "string", 'default': "multiple"},
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  #fastq文件夹
            {"name": "fq_type", "type": "string"},  #PE OR SE
            {"name": "group_table", "type": "infile", "format": "denovo_rna_v2.group_table"},
            {"name": "is_duplicate", "type": "bool", "default": ""},
            #有生物学重复的时候的分组文件
            {"name": "control_file", "type": "infile", "format": "denovo_rna_v2.compare_table"},
            #对照组文件，格式同分组文件
            {"name": "assembly_file", "type": "infile", 'format': "denovo_rna_v2.trinity_fasta"}, #trinity组装结果文件
            {"name": "gene_to_trans", "type": "infile", 'format': "denovo_rna_v2.common"}, #基因和转录本对应关系文件
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
        if not self.option('fastq_dir'):
            raise OptionError("必须输入fastq文件参数", code = "32003401")
        if not self.option('fastq_dir').prop['has_list_file']:
            raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件', code = "32003402")
        if self.option("sample_num") == "multiple":
            if not self.option('group_table'):
                raise OptionError("必须输入分组文件", code = "32003403")
            if not self.option('control_file'):
                raise OptionError("必须输入对照组文件", code = "32003404")
        if not self.option('fq_type').is_set:
            raise OptionError("必须设置测序类型：PE or SE", code = "32003405")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("测试类型只能是PE或者SE", code = "32003406")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "5G"


class FilecheckDenovoTool(Tool):
    """
    检查denovo rna输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(FilecheckDenovoTool, self).__init__(config)
        self.samples = list()

    def check_fastq(self):
        self.logger.info("正在检测fastq_dir文件")
        if not self.option("fastq_dir").prop["has_list_file"]:
            raise OptionError('fastq文件夹中必须含有一个名为list.txt的文件名--样本名的对应文件', code = "32003407")
        self.samples = self.option("fastq_dir").prop["samples"]
        for sample in self.samples:
            if len(sample) >= 20:
                raise OptionError("The length of sample name {} is longer than 20 characters.".format(sample))
            match = re.match(r'[^a-zA-Z]', sample)
            if match is not None:
                raise OptionError('Sample name {} must be start with english letter.'.format(sample))
            for char in sample:
                match = re.search(r'([^a-zA-z0-9_])', char)
                if match is not None:
                    raise OptionError('{} contains special character: {}'.format(sample, match.group()))
        ## 跳过fastq详细检查，减少流程运行时间
        '''
        col_num = self.get_list_info()
        if self.option("fq_type") in ["PE"] and col_num != 3:
            raise OptionError("PE文件夹的list应该包含三行信息，文件名-样本名-左端OR右端")
        file_list = FileSampleFile()
        list_txt = os.path.join(self.option('fastq_dir').prop['path'], "list.txt")
        file_list.set_path(list_txt)
        file_sample = file_list.get_list()
        # self.logger.info('%s' % file_sample)
        if self.option('fq_type') == 'PE':
            for i in file_sample.keys():
                if len(file_sample[i]) != 2:
                    raise OptionError("PE测序时，每个样本至少有一个左端fq和右端fq文件")
        files = self.option('fastq_dir').prop['fastq_basename']
        # self.logger.info('%s' % files)
        for f in files:
            fq_path = os.path.join(self.option('fastq_dir').prop['path'], f)
            my_fastq = FastqFile()
            my_fastq.set_path(fq_path)
            if re.search('\.gz$', f) or re.search('\.gzip$', f):
                my_fastq.is_gz = True
            my_fastq.check_content()
        '''
        self.logger.info("fastq文件检测完毕")

    def get_list_info(self):
        list_path = self.option("fastq_dir").prop["path"] + "/list.txt"
        with open(list_path, "r") as l:
            col_num = len(l.readline().strip().split())
        return col_num

    def check_group(self):
        if self.option('group_table').is_set:
            self.logger.info("正在检测group文件")
            self.option("group_table").get_info()
            gp_sample = self.option("group_table").prop["sample"]
            self.logger.info("group_table中有{}个样本".format(str(len(set(gp_sample)))))
            self.logger.info("fastq_dir中list有{}个样本".format(str(len(set(self.samples)))))
            if set(gp_sample) != set(self.samples):
                self.set_error("group和fastq_dir中list文件不对应")
            for gp in gp_sample:
                if gp not in self.samples:
                    self.set_error("group表出错, 样本%s在fastq文件中未出现", variables = (gp), code = "32003408")
        else:
            self.logger.info("未检测到group文件， 跳过...")
        if self.option("sample_num") == "multiple" and self.option("is_duplicate"):
            group_dict = self.option("group_table").prop['group_dict']
            group_names = set(group_dict.keys())
            samples = set(self.option("group_table").prop["sample"])
            coflict_names = group_names & samples
            if len(coflict_names) > 0 :
                self.set_error("group表出错, 当有生物学重复时,不允许组名和样本名相同")

        self.logger.info("group文件检测完毕")

    def check_control(self):
        self.logger.info("正在检测control文件")
        vs_list = self.option("control_file").prop["cmp_list"]
        con_samples = []
        if self.option('group_table').is_set:
            group_scheme = self.option('group_table').prop['group_scheme'][0]
            group_name = self.option('group_table').get_group_name(group_scheme)
        for vs in vs_list:
            for i in vs:
                if i not in con_samples:
                    con_samples.append(i)
            for cp in con_samples:
                if self.option('group_table').is_set:
                    if cp not in group_name:
                        self.set_error("control表出错，分组%s在fastq文件中未出现", variables = (cp), code = "32003409")
                else:
                    if cp not in self.samples:
                        self.set_error("control表出错，样本%s在fastq文件中未出现", variables = (cp), code = "32003410")
        self.logger.info("control文件检测完毕")

    def check_fasta(self):
        if self.option('assembly_file').is_set:
            self.logger.info("正在检测fasta文件")
            if self.option("assembly_file").prop["file_format"]!= 'FASTA':
                raise FileError("文件格式错误", code = "32003411")
            if self.option("assembly_file").prop["seq_number"] < 1:
                raise FileError("应该至少含有一条序列", code = "32003412")
        else:
            self.logger.info("未检测到fasta文件， 跳过...")
        self.logger.info("fasta文件检测完毕")


    def run(self):
        super(FilecheckDenovoTool, self).run()
        self.check_fastq()
        if self.option("sample_num") == "multiple":
            self.check_group()
            if self.option('control_file').is_set:
                self.check_control()
        self.check_fasta()
        self.end()
