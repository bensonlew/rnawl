# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/15'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir
import pandas as pd


class BacAssemble3Module(Module):
    """
    细菌基因组流程做二代数据拼接总模块
    """

    def __init__(self, work_id):
        super(BacAssemble3Module, self).__init__(work_id)
        option = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir", "required": True},
            {"name": "bam_dir", "type": "infile", "format": "bacgenome.simple_dir"},
            {"name": "assem_tool", "type": "string", "choose": ['canu', 'hgap/falcon'], "default": "canu"},
            {"name": "fq_info", "type": "infile", "format": "bacgenome.simple_file", "required": True},
            {"name": "bam_info", "type": "infile", "format": "bacgenome.simple_file"},
            {"name": "error_rate", "type": "string", "default": "0.013"},  # canu param
            {"name": "canu_version", "type": "string", "default": "1.3", 'choose': ["1.3", "1.7"]},  # canu version
            {"name": "corMinCoverage", "type": "string", "default": "default", 'choose': ['default', '0', '1', '2', "3", "4"]},  # canu param
            {"name": "corMhapSensitivity", "type": "string", "default": "default", "choose": ["default", "low", "normal", "high"]},  # canu param
            {"name": "output1", "type": "outfile", "format": "sequence.fasta_dir"},
            {"name": "output2", "type": "outfile", "format": "sequence.fasta_dir"}
            # canu的错误率参数
            # spades不在这里做，需要等抽取的pe reads完成
        ]
        self.add_option(option)
        self.first_assem_tools = []  # 首选的拼接工具，其结果用于最后的统计评估
        self.second_assem_tools = []  # 备选的拼接工具，结果用于交互分析
        self.sample_lib_dic = {}

    def check_options(self):
        """
        检查参数
        :return:
        """
        if self.option("bam_info").is_set:
            if not self.option("bam_dir").is_set:
                raise OptionError("输入bam_info必须提供bam_dir参数")
        data = pd.read_table(self.option("fq_info").prop["path"], index_col="Sample Name")
        self.sample_lib_dic = data["Library"].to_dict()
        return True


    def run(self):
        super(BacAssemble3Module, self).run()
        if self.option("assem_tool") == "canu":
            self.run_canu(first=True)
            self.run_hgap_falcon(first=False)
        elif self.option("assem_tool") == "hgap/falcon":
            self.run_hgap_falcon(first=True)
            self.run_canu(first=False)
        self.on_rely(self.first_assem_tools + self.second_assem_tools, self.set_output)
        for tool in self.first_assem_tools:
            tool.run()
        for tool in self.second_assem_tools:
            tool.run()

    def run_canu(self, first):
        data = pd.read_table(self.option("fq_info").prop["path"])
        data.apply(self.parse_canu, axis=1, first=first)

    def parse_canu(self, df, first):
        tool = self.add_tool("bacgenome.canu_assemble")
        opts = {
            "subread_fq": self.option("fq_dir").prop["path"] + "/" + df["File"],
            "genomeSize": str(df["genome"]),
            "sample_name": df["Sample Name"],
            "correctedErrorRate": self.option("error_rate"),
            "corMinCoverage": self.option("corMinCoverage"),
            "corMhapSensitivity": self.option("corMhapSensitivity"),
            "canu_version": "1.7"
        }
        if df["Library"] == "Pacbio":
            opts["read_type"] = "pacbio-raw"
        elif df["Library"] == "Nanopore":
            first = True  # Nanopore只能做canu拼接
            opts["read_type"] = "nanopore-raw"
        tool.set_options(opts)
        if first:
            self.first_assem_tools.append(tool)
        else:
            self.second_assem_tools.append(tool)

    def run_hgap_falcon(self, first):
        data = pd.read_table(self.option("fq_info").prop["path"])
        if self.option("bam_info").is_set:
            bam_data = pd.read_table(self.option("bam_info").prop["path"])
        else:
            bam_data = False
        data.apply(self.parse_hgap, axis=1, first=first, bam=bam_data)

    def parse_hgap(self, df, first, bam):
        sample_name = df["Sample Name"]
        # genome_size = df["genome"]  #  临时改变大小
        genome_size = 3
        file_name = df["File"]
        if isinstance(bam, pd.DataFrame):
            bam_samples = bam["Sample Name"].tolist()
            if sample_name in bam_samples:
                file_name = bam[bam["Sample Name"] == sample_name].iloc[0]["File"]
                self.add_hgap_tool(sample_name, genome_size, file_name, first)
            elif self.sample_lib_dic[sample_name] == "Pacbio":  # 即使只有fastq序列，nanopore数据运行不了falcon拼接
                self.add_falcon_tool(sample_name, genome_size, file_name, first)
        elif self.sample_lib_dic[sample_name] == "Pacbio":
            self.add_falcon_tool(sample_name, genome_size, file_name, first)

    def add_hgap_tool(self, sample_name, genome_size, file_name, first):
        tool = self.add_tool("assemble.hgap")
        tool.set_options({
            "bam_file": self.option("bam_dir").prop["path"] + "/" + file_name,
            "genome_size": "%s000000" % genome_size,
            "sample_name": sample_name
        })
        if first:
            self.first_assem_tools.append(tool)
        else:
            self.second_assem_tools.append(tool)

    def add_falcon_tool(self, sample_name, genome_size, file_name, first):
        tool = self.add_tool("assemble.falcon")
        tool.set_options({
            "fq_file": self.option("fq_dir").prop["path"] + "/" + file_name,
            "genome_size": "%s000000" % genome_size,
            "sample_name": sample_name
        })
        if first:
            self.first_assem_tools.append(tool)
        else:
            self.second_assem_tools.append(tool)

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for tool in self.first_assem_tools:
            link_dir(tool.output_dir, self.output_dir + "/first")
        for tool in self.second_assem_tools:
            link_dir(tool.output_dir, self.output_dir + "/second")
        try:
            self.option("output1").set_path(self.output_dir + "/first")
            self.option("output1").check()
        except:
            self.logger.info("first路径为空")
        try:
            self.option("output2").set_path(self.output_dir + "/second")
            self.option("output2").check()
        except:
            self.logger.info("second路径为空")
        self.logger.info(self.option("output1").is_set)
        self.logger.info(self.option("output2").is_set)
        self.logger.info("设置组装结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BacAssemble3Module, self).end()
