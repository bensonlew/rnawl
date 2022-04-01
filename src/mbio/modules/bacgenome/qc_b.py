# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/9'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd
from mbio.packages.metagenomic.common import link_dir,link_file


class QcBModule(Module):
    """
    fastp的质控module
    """

    def __init__(self, work_id):
        super(QcBModule, self).__init__(work_id)
        option = [
            {"name": "list", "type": "infile", "format": "bacgenome.simple_file"},
            {"name": "fastq_dir", "type": "infile", "format": "bacgenome.simple_dir"},
            {"name": "phix_tool", "type": "string", "default": "bwa", "choose": ["bwa", "bowtie"]}
        ]
        self.add_option(option)
        self.raw_stat_tool  = self.add_module("bacgenome.hiseq_seq_stat")
        self.qc_tool = self.add_module("bacgenome.bac_high_seq_qc")
        self.qc_stat_tool = self.add_module("bacgenome.hiseq_seq_stat")
        self.fastq_dir = os.path.join(self.work_dir, "qc_input")
        self.stat_info = pd.DataFrame(columns=["Sample_lib", "pair reads(#)", "total bases(bp)", "Q20(%)", "Q30(%)", "type"])

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(QcBModule, self).run()
        link_dir(self.option("fastq_dir").prop["path"], self.fastq_dir)
        self.on_rely([self.raw_stat_tool, self.qc_tool], self.run_qc_stat)
        self.qc_stat_tool.on("end", self.set_output)
        self.run_raw_stat()
        self.run_qc()



    def get_list(self):
        self.raw_list = os.path.join(self.fastq_dir, "list.txt")
        self.sample_info = os.path.join(self.fastq_dir, "sample_info")
        data = pd.read_table(self.option("list").prop["path"], header=None)
        # data["sample_lib"] = data[0] + "_" + data[1] + "_" + data[2]
        data["sample_lib"] = data.apply(lambda x: x[6][:-5], axis=1)
        data[[6, "sample_lib", 5]].to_csv(self.raw_list, sep="\t", header=False, index=False)
        data[["sample_lib", 3, 4]].to_csv(self.sample_info, sep="\t", header=False, index=False)

    def run_raw_stat(self):
        self.get_list()
        opts = {
            "list": self.raw_list,
            "fastq_dir": self.option("fastq_dir")
        }
        self.raw_stat_tool.set_options(opts)
        self.raw_stat_tool.run()

    def run_qc(self):
        self.qc_tool.set_options({
            'sample_path': self.raw_list,
            'sample_info': self.sample_info,
            'rm_single': True,
            'phix_tool': self.option("phix_tool")
        })
        self.qc_tool.run()

    def run_qc_stat(self):
        self.qc_stat_tool.set_options({
            "list": self.qc_tool.option("clean_list").path,
            "fastq_dir": self.qc_tool.output_dir
        })
        self.qc_stat_tool.run()

    def parse_stat(self, file_path, type):
        data = pd.read_table(file_path)
        data["type"] = type
        self.stat_info = self.stat_info.append(data)

    def write_stat(self):
        stat_path = os.path.join(self.output_dir, "stat")
        if not os.path.isdir(stat_path):
            os.mkdir(stat_path)
        clean_data = self.stat_info[self.stat_info["type"]=="clean"].drop("type", axis=1)
        raw_data = self.stat_info[self.stat_info["type"]=="raw"].drop("type",axis=1)
        clean_data.to_csv(stat_path + "/clean_statistics.xls", sep="\t", header=True, index=False)
        raw_data.to_csv(stat_path + "/raw_statistics.xls", sep="\t", header=True, index=False)

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        fastx_path = os.path.join(self.output_dir, "fastx")
        data_path = os.path.join(self.output_dir, "cleandata")
        if not os.path.isdir(fastx_path):
            os.mkdir(fastx_path)
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        for tool in self.raw_stat_tool.stat_tools:
            files = os.listdir(tool.output_dir)
            for file in files:
                file_path = os.path.join(tool.output_dir, file)
                if file == "qc_stat.xls":
                    self.parse_stat(file_path, type="raw")
                else:
                    # "eg. normal_gz_PE_430.2.fq_fastxstat"
                    file_prefix = file[:-15]
                    if file[-14] == "1":
                        new_file = file_prefix + "_l.raw_fastxstat"
                    else:
                        new_file = file_prefix + "_r.raw_fastxstat"
                    new_file_path = os.path.join(fastx_path, new_file)
                    link_file(file_path, new_file_path)
        for tool in self.qc_stat_tool.stat_tools:
            files = os.listdir(tool.output_dir)
            for file in files:
                file_path = os.path.join(tool.output_dir, file)
                if file == "qc_stat.xls":
                    self.parse_stat(file_path, type="clean")
                else:
                    # "eg. normal_gz_PE_430.clean.2.fq_fastxstat"
                    file_prefix = file[:-21]
                    if file[-14] == "1":
                        new_file = file_prefix + "_l.clean_fastxstat"
                    else:
                        new_file = file_prefix + "_r.clean_fastxstat"
                    new_file_path = os.path.join(fastx_path, new_file)
                    link_file(file_path, new_file_path)
        link_dir(self.qc_tool.output_dir, data_path)
        self.logger.info(self.stat_info)
        self.write_stat()  # q20 q30 文本，需要上一级的module添加样本实际名称，insert_size, reads length等信息，并整合成统一文档
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
        super(QcBModule, self).end()
