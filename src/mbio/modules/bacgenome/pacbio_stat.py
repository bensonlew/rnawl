# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/10'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir,link_file
import pandas as pd


class PacbioStatModule(Module):
    """
    三代数据的序列长度统计程序
    """

    def __init__(self, work_id):
        super(PacbioStatModule, self).__init__(work_id)
        option = [
            {"name": "third_list", "type": "infile", "format": "bacgenome.simple_file", "required": False},
            {"name": "fastq_dir", "type": "infile", "format": "bacgenome.simple_dir", "required": True}
        ]
        self.add_option(option)
        self.run_tools = []
        self.stat_info = pd.DataFrame(columns=["Sample", "Total Reads No.", "Total Bases (bp)", "Largest (bp)", "Average Len (bp)"])

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(PacbioStatModule, self).run()
        self.run_stat()

    def run_stat(self):
        with open(self.option("third_list").prop["path"], "r") as file:
            for line in file.readlines():
                line = line.strip().split("\t")
                sample_name = line[0]
                input_fq = os.path.join(self.option("fastq_dir").prop["path"], line[2])
                tool = self.add_tool('fungi_genome.pacbio_clean')
                tool.set_options({
                    "input_fq": input_fq,
                    "sample_name": sample_name
                })
                self.run_tools.append(tool)
        if len(self.run_tools) == 1:
            self.run_tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.run_tools, self.set_output)
        for tool in self.run_tools:
            tool.run()

    def parse_stat(self, name, file_path):
        data = pd.read_table(file_path)
        data["Sample"] = name
        self.stat_info = self.stat_info.append(data)

    def write_stat(self, output):
        self.stat_info.reindex(columns=["Sample", "Total Reads No.", "Total Bases (bp)", "Largest (bp)", "Average Len (bp)"])\
            .to_csv(output, sep="\t", header=True, index=False)

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for tool in self.run_tools:
            name = tool.option("sample_name")
            for file in os.listdir(tool.output_dir):
                file_path = os.path.join(tool.output_dir, file)
                if file == name + ".clean.len.xls":
                    link_file(file_path, os.path.join(self.output_dir, name + ".len.xls"))
                elif file == name + ".PacBio_statistics.xls":
                    self.parse_stat(name, file_path)
        self.write_stat(os.path.join(self.output_dir, "statistics.xls"))  # 需要上级module加文库信息
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
        super(PacbioStatModule, self).end()
