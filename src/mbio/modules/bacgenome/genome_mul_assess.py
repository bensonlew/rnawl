# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/16'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd
from mbio.packages.metagenomic.common import link_dir,link_file


class GenomeMulAssessModule(Module):
    """
    多样本的基因组评估
    """

    def __init__(self, work_id):
        super(GenomeMulAssessModule, self).__init__(work_id)
        option = [
            {"name": "fastq_list", "type": "infile", "format": "bacgenome.simple_file"},
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir", "required": True},
            {"name": "scaf_dir", "type": "infile", "format": "sequence.fasta_dir", "required": True},
            {"name": "pe_sample_list", "type": "string"}, # 哪些样本是做二代数据的基因组评估流程
            {"name": "chr_sample_list", "type": "string"} # 只做pca的样品
        ]
        self.add_option(option)
        self.run_tools = []  # edit tool list
        self.samples = []
        self.chr_samples = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("pe_sample_list") and not self.option("chr_sample_list"):
            raise OptionError("pe_sample_list or chr_sample_list must be setted at least one")
        if self.option("pe_sample_list"):
            self.samples = self.option("pe_sample_list").split(",")
        if self.samples and not self.option("fastq_list").is_set:
            raise OptionError("必须输入fastq_list")
        if self.option("chr_sample_list"):
            self.chr_samples = self.option("chr_sample_list").split(",")
        return True

    def run(self):
        super(GenomeMulAssessModule, self).run()
        if not os.path.isdir(os.path.join(self.work_dir, "assess_input")):
            os.mkdir(os.path.join(self.work_dir, "assess_input"))
        if self.samples:
            self.run_assess()
        if self.chr_samples:
            self.run_simple_assess()
        self.on_rely(self.run_tools, self.set_output)
        for tool in self.run_tools:
            tool.run()

    def run_assess(self):
        data = pd.read_table(self.option("fastq_dir").prop["path"] + "/list.txt", header=None)
        with open(self.option("fastq_list").prop["path"], "r") as file:
            lines = file.readlines()[1:]
            for line in lines:
                sample_list, base_num, genome_size = line.strip().split("\t")
                sample = sample_list.split("_PE_")[0]
                if sample not in self.samples:
                    continue
                scaf = sample + "_scaf.fna"
                assess_tool = self.add_module("bacgenome.genome_assess2")
                self.logger.info(os.path.join(self.option("scaf_dir").prop["path"], scaf))
                assess_tool.set_options({
                    "seq": os.path.join(self.option("scaf_dir").prop["path"], scaf),
                    "fastq_dir": self.get_fastq(sample_list, data),
                    "bases": base_num,
                    "sample_name": sample
                })
                self.run_tools.append(assess_tool)

    def run_simple_assess(self):
        # 因为没有fastq，所以不做gc_depth等分析，只做pca和数据库比对
        for sample in self.chr_samples:
            scaf = sample + "_scaf.fna"
            assess_tool = self.add_module("bacgenome.genome_assess2")
            assess_tool.set_options({
                "seq": os.path.join(self.option("scaf_dir").prop["path"], scaf),
                "sample_name": sample
            })
            self.run_tools.append(assess_tool)

    def get_fastq(self, sample_list, data):
        sample = sample_list.split("_PE_")[0]
        sample_fq_path = os.path.join(self.work_dir, "assess_input", sample)
        if not os.path.isdir(sample_fq_path):
            os.mkdir(sample_fq_path)
        fq_info = data[data[1]==sample_list]
        fq_info.iloc[:,1] = sample
        fq1 = data[(data[1]==sample_list) & (data[2]=="l")].iloc[0,0]
        fq2 = data[(data[1]==sample_list) & (data[2]=="r")].iloc[0,0]
        link_file(os.path.join(self.option("fastq_dir").prop["path"], fq1), os.path.join(sample_fq_path, fq1))
        link_file(os.path.join(self.option("fastq_dir").prop["path"], fq2), os.path.join(sample_fq_path, fq2))
        fq_info.to_csv(os.path.join(sample_fq_path, "list.txt"), sep="\t", header=False, index=False)
        return sample_fq_path

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for tool in self.run_tools:
            sample_name = tool.option("sample_name")
            link_dir(tool.output_dir, self.output_dir + "/" + sample_name)
        self.logger.info("设置注释结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(GenomeMulAssessModule, self).end()
