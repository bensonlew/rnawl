# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/5/6'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd
from Bio import SeqIO
from mbio.packages.metagenomic.common import link_file

class ChrGapStepModule(Module):
    """
    做一种数据的完成图结果gap矫正
    """

    def __init__(self, work_id):
        super(ChrGapStepModule, self).__init__(work_id)
        option = [
            {"name": "chr_fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "fill_fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "pre_log", "type": "infile", "format": "sequence.profile_table"},
            {"name": "overlap", "type": "int", "default": 1000},
            {"name": "identity", "type": "float", "default": 0.95},
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "log", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(option)
        self.run_tools = []  # edit tool list
        self.all_circled = False
        self.log_have_circle_info = False
        self.split_fa = self.add_tool('bacgenome.split_chr')  #拆分骨架序列环性和线性分开
        self.gap_fa = self.add_tool("bacgenome.gap_fa")  # 补gap

    def check_options(self):
        """
        检查参数
        :return:
        """
        # edit options check
        return True

    def run(self):
        super(ChrGapStepModule, self).run()
        self.split_fa.on("end", self.run_gap)
        self.gap_fa.on("end", self.set_output)
        self.run_split_chr()

    def run_split_chr(self):
        self.split_fa.set_options({
            "fa": self.option("chr_fa"),
            "pre_log": self.option("pre_log")
        })
        self.split_fa.run()

    def run_gap(self):
        pre_log = self.option("pre_log")
        fasta = self.get_fasta(self.option("fill_fa").prop['path'])
        opts = {
            "chr_fa": self.split_fa.option("line"),
            "fill_fa": fasta,
            "pre_log": pre_log,
            "overlap": self.option("overlap"),
            "identity": self.option("identity")
        }
        if self.split_fa.option("circle").is_set:
            opts["circle_fa"] = self.split_fa.option("circle")
        self.gap_fa.set_options(opts)
        self.gap_fa.run()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        new_log = os.path.join(self.output_dir, "result.log")
        data = pd.read_table(self.gap_fa.option("log").prop["path"], header=None)
        if data.loc[~data[3],].empty:
            self.all_circled = True
        seq_records = SeqIO.parse(self.gap_fa.option("out_fa").prop["path"], "fasta")
        for record in seq_records:
            data.loc[data[1]==record.id, 4] = str(len(record.seq))
        data[3] = data[3].map({True: "Circular", False: "Linear"})
        data.to_csv(new_log, sep="\t", header=False, index=False)
        self.option("out_fa").set_path(self.gap_fa.option("out_fa").prop["path"])
        self.option("log").set_path(new_log)
        link_file(self.gap_fa.output_dir + "/blast.xls", self.output_dir + "/blast.xls")
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
        super(ChrGapStepModule, self).end()

    def get_fasta(self, fa):
        list =[]
        for seq_record in SeqIO.parse(fa, "fasta"):
            if len(seq_record) >= 2000:
                list.append(seq_record)
        SeqIO.write(list, self.work_dir + "/seq2.fasta", "fasta")
        return self.work_dir + "/seq2.fasta"