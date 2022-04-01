# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/15'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_dir
import pandas as pd


class BacAssemble2Module(Module):
    """
    细菌基因组流程做二代数据拼接总模块
    """

    def __init__(self, work_id):
        super(BacAssemble2Module, self).__init__(work_id)
        option = [
            {"name": "fq_dir", "type": "infile", "format": "sequence.fastq_dir", "required": True},
            {"name": "assem_tool", "type": "string", "choose": ['soapdenovo', 'velvet', 'spades'], "default": "soapdenovo"},
            {"name": "major_info", "type": "infile", "format": "bacgenome.simple_file", "required": True},
            {"name": "sample_info", "type": "infile", "format": "bacgenome.simple_file", "required": True},
            {"name": "kmers", "type": "string"}
            # 其他的拼接参数一大堆......
        ]
        soapdenovo_opt = [
            {"name": "soapdenovo_D", "type": "int", "default": 1, "min": 1, "max": 10},
            {"name": "soapdenovo_d", "type": "string", "default": "3,5,10"},
            {"name": "soapdenovo_M", "type": "int", "default": 1, "min":0, "max": 3},
            {"name": "soapdenovo_R", "type": "bool", "default": True},
            {"name": "soapdenovo_F", "type": "bool", "default": True},
            {"name": "soapdenovo_u", "type": "string", "default": "unmask", "choose": ["mask", "unmask"]},
            {"name": "soapdenovo_G", "type": "int", "default": 50, "min": 0}
        ]
        velvet_opt = [
            {"name": "velvet_min_contig_lgth", "type": "int", "default": 200, "min": 100},
            {"name": "velvet_min_pair_count", "type": "int", "default": 15, "min": 5}
        ]
        option += soapdenovo_opt + velvet_opt
        self.add_option(option)
        self.assem_tools = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        return True

    def run(self):
        super(BacAssemble2Module, self).run()
        self.run_assemble()

    def run_assemble(self):
        with open(self.option("major_info").prop["path"], "r") as file:
            lines = file.readlines()[1:]
            for line in lines:
                line = line.strip().split("\t")
                sample_name = line[0].split("_PE_")[0]
                assemble = self.add_module("bacgenome.assemble_hiseq")
                assemble.set_options({
                    "fq_dir": self.option("fq_dir"),
                    "sample_name": sample_name,
                    "assem_tool": self.option("assem_tool"),
                    "sample_info": self.option("sample_info"),
                    "kmers": self.option("kmers"),
                    "soapdenovo_D": self.option("soapdenovo_D"),
                    "soapdenovo_d": self.option("soapdenovo_d"),
                    "soapdenovo_M": self.option("soapdenovo_M"),
                    "soapdenovo_R": self.option("soapdenovo_R"),
                    "soapdenovo_F": self.option("soapdenovo_F"),
                    "soapdenovo_u": self.option("soapdenovo_u"),
                    "soapdenovo_G": self.option("soapdenovo_G"),
                    "velvet_min_contig_lgth": self.option("velvet_min_contig_lgth"),
                    "velvet_min_pair_count": self.option("velvet_min_pair_count")
                })
                self.assem_tools.append(assemble)
        self.on_rely(self.assem_tools, self.set_output)
        for tool in self.assem_tools:
            tool.run()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for tool in self.assem_tools:
            link_dir(tool.output_dir, self.output_dir)
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
        super(BacAssemble2Module, self).end()
