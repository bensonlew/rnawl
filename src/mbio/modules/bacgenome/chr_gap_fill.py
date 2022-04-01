# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/5/6'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagenomic.common import link_file


class ChrGapFillModule(Module):
    """
    交互分析校正序列模块
    """

    def __init__(self, work_id):
        super(ChrGapFillModule, self).__init__(work_id)
        option = [
            {"name": "chr_fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "cir_fa", "type": "infile", "format": "sequence.fasta"},
            {"name": "pre_log", "type": "infile", "format": "sequence.profile_table"},  # 前一分析记录中的log，对于origin的记录，有chr_fa结果用到的拼接软件
            {"name": "genome_json", "type": "string", "required": True},
            {"name": "overlap", "type": "int", "default": 1000},
            {"name": "identity", "type": "float", "default": 0.95},
            {"name": "log", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"}
        ]
        self.add_option(option)
        self.run_tools = []  # edit tool list
        self.fasta_needed = self.add_tool("bacgenome.get_fa_from_json")
        self.gap = self.add_module("bacgenome.chr_gap_step")
        self.gap_data_list = []  # pe ,falcon/hgap, spades, canu

    def check_options(self):
        """
        检查参数
        :return:
        """
        # edit options check
        if not self.option("cir_fa").is_set:
            raise OptionError("沒有可用于校正的数据")
        return True

    def run(self):
        super(ChrGapFillModule, self).run()
        self.get_fa_from_json()

    def get_fa_from_json(self):
        self.fasta_needed.set_options({
            "chr_fa": self.option("chr_fa"),
            "genome_json": self.option("genome_json"),
            "log": self.option("pre_log")
        })
        self.fasta_needed.on("end", self.run_gap)
        self.fasta_needed.run()

    def run_gap(self):
        self.gap_data_list.append("gap assemble result")
        self.gap.set_options({
            "chr_fa": self.fasta_needed.option("out_fa"),# self.option("chr_fa"),
            "fill_fa": self.option("cir_fa"),
            "pre_log": self.fasta_needed.option("out_log"),
            "overlap": self.option("overlap"),
            "identity": self.option("identity")
        })  # module结果中包含一个fasta,和对应的log(每条序列来自原始的哪一条序列)
        # module包含all_circled变量,默认为False
        self.gap.on("end", self.set_output)
        self.gap.run()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        link_file(self.gap.option("out_fa").prop["path"], self.output_dir + "/result.fa")
        link_file(self.gap.option("log").prop["path"], self.output_dir + "/result.log")
        link_file(self.gap.output_dir + "/blast.xls", self.output_dir + "/blast.xls")
        self.option("out_fa").set_path(self.output_dir + "/result.fa")
        self.option("log").set_path(self.output_dir + "/result.log")
        self.logger.info("设置结果成功")
        self.end()

    def end(self):
        super(ChrGapFillModule, self).end()
