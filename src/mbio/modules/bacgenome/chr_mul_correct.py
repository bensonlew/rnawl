# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/28'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd
from mbio.packages.metagenomic.common import link_file


class ChrMulCorrectModule(Module):
    """
    对完成图有pe reads的结果做pe_reads校正和spades拼接，
    全部去做，无论是否完全成环，可能里面有gap区，还是要补的
    """

    def __init__(self, work_id):
        super(ChrMulCorrectModule, self).__init__(work_id)
        option = [
            {"name": "pe_dir", "type": "infile", "format": "sequence.fastq_dir", "required": True},  # 用质控后的结果
            {"name": "chr_dir", "type": "infile", "format": "sequence.fastq_dir", "required": True},  # 原始数据解压路径
            {"name": "major_info", "type": "infile", "format": "bacgenome.simple_file", "required": True},  # 最多base数文库的pe列表
            {"name": "seq_scaf1", "type": "infile", "format": "sequence.fasta_dir", "required": True},  # 首选的三代结果
            {"name": "chr_info", "type": "infile", "format": "bacgenome.simple_file", "required": True}, # 三代fq数据的信息
            {"name": "seq_scaf2", "type": "infile", "format": "sequence.fasta_dir", "required": True}, # 备选的三代结果
            {"name": "pe_sample_list", "type": "string", "required": True}
        ]
        self.add_option(option)
        self.run_tools = []
        self.spades = []
        self.have_pe_samples = []
        self.pe_list = {}
        self.chr_path = {}
        self.lib_type = {}
        self.fq1_path = {}
        self.fq2_path = {}

    def check_options(self):
        """
        检查参数
        :return:
        """
        self.have_pe_samples = self.option("pe_sample_list").split(",")
        return True

    def run(self):
        super(ChrMulCorrectModule, self).run()
        self.get_info()
        self.run_spades()
        self.run_correct_major()
        self.run_correct_minor()
        self.on_rely(self.run_tools + self.spades, self.set_output)
        for tool in self.run_tools:
            tool.run()
        for tool in self.spades:
            tool.run()

    def get_info(self):
        major_data = pd.read_table(self.option("major_info").prop["path"])
        major_data.apply(self.filter_pe_sample, axis=1)
        chr_data = pd.read_table(self.option("chr_info").prop["path"])
        chr_data.apply(self.get_chr, axis=1)

    def filter_pe_sample(self, df):
        sample_lib = df["Sample_lib"]
        sample_name = sample_lib.split("_PE_")[0]
        if sample_name not in self.have_pe_samples:
            return 1
        self.fq1_path[sample_name] = os.path.join(self.option("pe_dir").prop["path"], sample_lib + ".clean.1.fq")
        self.fq2_path[sample_name] = os.path.join(self.option("pe_dir").prop["path"], sample_lib + ".clean.2.fq")
        self.write_pe_list(sample_name)
        return 1

    def write_pe_list(self, sample_name):
        self.pe_list[sample_name] = self.work_dir + "/%s_list" % sample_name
        with open(self.pe_list[sample_name], "w") as file:
            file.write("%s\t%s\t%s\t\t1\n" % (sample_name, self.fq1_path[sample_name], self.fq2_path[sample_name]))

    def get_chr(self, df):
        self.lib_type[df["Sample Name"]] = df["Library"]
        self.chr_path[df["Sample Name"]] = os.path.join(self.option("chr_dir").prop["path"], df["File"])
        return 1

    def run_spades(self):
        for sample in self.have_pe_samples:
            tool = self.add_tool("assemble.assemble_spades")
            opts = {
                "PE_list" : self.pe_list[sample],
                "sample_name": sample,
            }
            if self.lib_type[sample] == "Nanopore":
                opts["nanopore"] = self.chr_path[sample]
            else:
                opts["pacbio"] = self.chr_path[sample]
            tool.set_options(opts)
            self.spades.append(tool)

    def run_correct_major(self):
        for file in os.listdir(self.option("seq_scaf1").prop["path"]):
            full_path = os.path.join(self.option("seq_scaf1").prop["path"], file)
            sample = file.split(".scaf")[0]
            if sample not in self.have_pe_samples:
                continue
            tool = self.add_module("bacgenome.pilon_correct")
            tool.set_options({
                "seq_scaf": full_path,
                "fq1": self.fq1_path[sample],
                "fq2": self.fq2_path[sample],
                "is_major_result": True,
                "sample_name": sample
            })
            self.run_tools.append(tool)

    def run_correct_minor(self):
        for file in os.listdir(self.option("seq_scaf2").prop["path"]):
            full_path = os.path.join(self.option("seq_scaf2").prop["path"], file)
            sample = file.split(".scaf")[0]
            if sample not in self.have_pe_samples:
                continue
            tool = self.add_module("bacgenome.pilon_correct")
            tool.set_options({
                "seq_scaf": full_path,
                "fq1": self.fq1_path[sample],
                "fq2": self.fq2_path[sample],
                "is_major_result": False,
                "sample_name": sample
            })
            self.run_tools.append(tool)

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        first_dir = os.path.join(self.output_dir, "first")
        second_dir = os.path.join(self.output_dir, "second")
        third_dir = os.path.join(self.output_dir, "third")
        coverage_dir = os.path.join(self.output_dir, "coverage")
        for dir in [first_dir, second_dir, third_dir, coverage_dir]:
            if not os.path.isdir(dir):
                os.mkdir(dir)
        for tool in self.spades:
            file_name = os.path.basename(tool.option("scf_seq").prop["path"])
            link_file(tool.option("scf_seq").prop["path"], os.path.join(third_dir, file_name))
        for tool in self.run_tools:
            if tool.option("is_major_result"):
                link_file(tool.option("scaffold").prop["path"], os.path.join(first_dir, tool.option("sample_name") + ".scaffold.fna"))
                link_file(tool.option("coverage").prop["path"], os.path.join(coverage_dir, tool.option("sample_name") + ".abund"))
                link_file(tool.option("mean_cov").prop["path"], os.path.join(coverage_dir, tool.option("sample_name") + ".depth"))
            else:
                link_file(tool.option("scaffold").prop["path"], os.path.join(second_dir, tool.option("sample_name") + ".scaffold.fna"))
        self.logger.info("设置结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(ChrMulCorrectModule, self).end()
