# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/26'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class PilonCorrectModule(Module):
    """
    
    """

    def __init__(self, work_id):
        super(PilonCorrectModule, self).__init__(work_id)
        options = [
            {"name": "seq_scaf", "type": "infile", "format": "sequence.fasta"},  # 最佳Gapcloser的scaffold文件
            # {'name': 'sample_name', "type": "string"},  # 样本名
            {"name": "is_major_result", "type": "bool", "required": True}, # 是否是主要的三代拼接结果，如果是，则计算coverage
            {'name': "fq1", "type": "infile", "format": "sequence.fastq", "required": True},
            {"name": "fq2", "type": "infile", "format": "sequence.fastq", "required": True},
            {"name": "sample_name", "type": "string"},
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
            {"name": "coverage", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "mean_cov", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self.correct_one = self.add_module("bacgenome.mapping")
        self.correct_two = self.add_module("bacgenome.mapping")
        self.correct_three = self.add_module("bacgenome.mapping")
        self.correct_four = self.add_module("bacgenome.mapping")
        self.correct_five = self.add_module("bacgenome.mapping")
        self.modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        # edit options check
        return True

    def run_correct_one(self):
        self.correct_one.set_options({
            "ref_fa": self.option("seq_scaf"),
            "fq1": self.option("fq1"),
            "fq2": self.option("fq2"),
            "correct_reads": True,
            "coverage": self.option("is_major_result"),
            "cal_cov_auto": self.option("is_major_result")
        })
        self.correct_one.on("end", self.run_correct_two)
        self.modules.append(self.correct_one)
        self.correct_one.run()

    def run_correct_two(self):
        if self.correct_one.correct_reads_is_complete:
            self.set_output()
            return
        self.correct_two.set_options({
            "ref_fa": self.correct_one.option("scf_seq"),
            "fq1": self.option("fq1"),
            "fq2": self.option("fq2"),
            "correct_reads": True,
            "coverage": self.option("is_major_result"),
            "cal_cov_auto": self.option("is_major_result")
        })
        self.correct_two.on("end", self.run_correct_three)
        self.modules.append(self.correct_two)
        self.correct_two.run()

    def run_correct_three(self):
        if self.correct_two.correct_reads_is_complete:
            self.set_output()
            return
        self.correct_three.set_options({
            "ref_fa": self.correct_two.option("scf_seq"),
            "fq1": self.option("fq1"),
            "fq2": self.option("fq2"),
            "correct_reads": True,
            "coverage": self.option("is_major_result"),
            "cal_cov_auto": self.option("is_major_result")
        })
        self.correct_three.on("end", self.run_correct_four)
        self.modules.append(self.correct_three)
        self.correct_three.run()

    def run_correct_four(self):
        if self.correct_three.correct_reads_is_complete:
            self.set_output()
            return
        self.correct_four.set_options({
            "ref_fa": self.correct_three.option("scf_seq"),
            "fq1": self.option("fq1"),
            "fq2": self.option("fq2"),
            "correct_reads": True,
            "coverage": self.option("is_major_result"),
            "cal_cov_auto": self.option("is_major_result")
        })
        self.correct_four.on("end", self.run_correct_five)
        self.modules.append(self.correct_four)
        self.correct_four.run()

    def run_correct_five(self):
        if self.correct_four.correct_reads_is_complete:
            self.set_output()
            return
        self.correct_five.set_options({
            "ref_fa": self.correct_four.option("scf_seq"),
            "fq1": self.option("fq1"),
            "fq2": self.option("fq2"),
            "correct_reads": True,
            "coverage": self.option("is_major_result"),
            "cal_cov_auto": False
        })
        self.correct_five.on("end", self.set_output)
        self.modules.append(self.correct_five)
        self.correct_five.run()

    def run(self):
        super(PilonCorrectModule, self).run()
        self.run_correct_one()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        last_result = ""
        for tool in self.modules:
            if tool.correct_reads_is_complete:
                last_result = tool
                break
        last_result = last_result if last_result else self.correct_five
        self.option("scaffold").set_path(last_result.option("scf_seq").prop["path"])
        if self.option("is_major_result"):
            self.option("coverage").set_path(last_result.option("cov_out").prop["path"])
            self.option("mean_cov").set_path(last_result.option("cov_out2").prop["path"])
        self.logger.info("设置pilon_correct结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(PilonCorrectModule, self).end()
