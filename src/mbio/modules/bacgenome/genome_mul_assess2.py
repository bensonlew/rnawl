# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# __modify__ = '2020/3/30'

import os,re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd
from mbio.packages.metagenomic.common import link_dir,link_file


class GenomeMulAssess2Module(Module):
    """
    多样本的基因组评估
    """

    def __init__(self, work_id):
        super(GenomeMulAssess2Module, self).__init__(work_id)
        option = [
            {"name": "fastq_list", "type": "infile", "format": "bacgenome.simple_file"},
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir", "required": True},
            {"name": "dir", "type": "string"},
            {"name": "pe_sample_list", "type": "string"}, # 哪些样本是做二代数据的基因组评估流程
        ]
        self.add_option(option)
        self.run_tools = []  # edit tool list
        self.samples = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("pe_sample_list"):
            raise OptionError("pe_sample_list must be setted at least one")
        if re.search(",", self.option("pe_sample_list")):
            self.samples = self.option("pe_sample_list").split(",")
        else:
            self.samples.append(self.option("pe_sample_list"))
        if self.samples and not self.option("fastq_list").is_set:
            raise OptionError("必须输入fastq_list")
        return True

    def run(self):
        super(GenomeMulAssess2Module, self).run()
        if not os.path.isdir(os.path.join(self.work_dir, "assess_input")):
            os.mkdir(os.path.join(self.work_dir, "assess_input"))
        self.run_assess()

    def run_assess(self):
        self.dict = {}
        self.dict_base, self.fq_dict = self.get_info(self.option('fastq_list').prop['path'])
        self.logger.info(self.samples)
        self.logger.info(self.dict_base)
        self.logger.info(self.fq_dict)
        for sample in self.samples:
            assess_tool = self.add_module("bacgenome.genome_assess3")
            self.dict[assess_tool] = sample
            assess_tool.set_options({
                "seq": self.option("dir") + "/" + sample + "/soapdenovo/" + sample + "_scaf.fna",
                "fastq_dir": self.fq_dict[sample],
                "bases": self.dict_base[sample],
                "sample_name": sample
            })
            self.run_tools.append(assess_tool)
        if len(self.run_tools) > 1:
            self.on_rely(self.run_tools, self.set_output)
        elif len(self.run_tools) == 1:
            self.run_tools[0].on('end', self.set_output)
        for module in self.run_tools:
            module.run()

    def get_info(self, input):
        dict = {}
        dict2 = {}
        data = pd.read_table(self.option("fastq_dir").prop["path"] + "/list.txt", header=None)
        with open(input, "r") as file:
            lines = file.readlines()[1:]
            for line in lines:
                sample_list, base_num, genome_size = line.strip().split("\t")
                sample = sample_list.split("_PE_")[0]
                dict[sample] = base_num
                sample_fq_path = os.path.join(self.work_dir, "assess_input", sample)
                if not os.path.isdir(sample_fq_path):
                    os.mkdir(sample_fq_path)
                fq_info = data[data[1] == sample_list]
                fq_info.iloc[:, 1] = sample
                fq1 = data[(data[1] == sample_list) & (data[2] == "l")].iloc[0, 0]
                fq2 = data[(data[1] == sample_list) & (data[2] == "r")].iloc[0, 0]
                link_file(os.path.join(self.option("fastq_dir").prop["path"], fq1), os.path.join(sample_fq_path, fq1))
                link_file(os.path.join(self.option("fastq_dir").prop["path"], fq2), os.path.join(sample_fq_path, fq2))
                fq_info.to_csv(os.path.join(sample_fq_path, "list.txt"), sep="\t", header=False, index=False)
                dict2[sample] = sample_fq_path
        return dict, dict2

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if len(self.run_tools) > 1:
            for tool in self.run_tools:
                link_dir(tool.output_dir, self.output_dir + "/" + self.dict[tool])
        elif len(self.run_tools) == 1:
            link_dir(self.run_tools[0].output_dir, self.output_dir + "/" + self.dict[self.run_tools[0]])
        self.logger.info("设置评估结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(GenomeMulAssess2Module, self).end()
