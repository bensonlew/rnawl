# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __modify__ = '2019/4/22'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import pandas as pd


class MappingModule(Module):
    """
    做pe reads mapping，根据需要计算coverage
    """

    def __init__(self, work_id):
        super(MappingModule, self).__init__(work_id)
        option = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta", "required": True},  # 组装结果
            {"name": "fq1", "type": "infile", "format": "sequence.fastq", "required": True},  # pe left reads
            {"name": "fq2", "type": "infile", "format": "sequence.fastq", "required": True},  # pe right reads
            {"name": "correct_reads", "type": "bool", "default": False},  # 是否根据比对的结果做ref_fa的校正
            {"name": "coverage", "type": "bool", "default": True},  # 是否计算coverage
            {"name": "cal_cov_auto", "type": "bool", "default": False},  # 是否自动判断是否计算coverage,会覆盖coverage参数
            {"name": "bam_out", "type": "outfile", "format": "align.bwa.bam"},
            {"name": "sam_out", "type": "outfile", "format": "align.bwa.sam"},
            {"name": "cov_out", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "cov_out2", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "scf_seq", "type": "outfile", "format": "sequence.fasta"}
        ]
        self.add_option(option)
        self.run_tools = []  # edit tool list
        self.bowtie2 = self.add_tool("align.bowtie2")
        self.convert = self.add_tool('metagbin.convert_format')
        self.coverage = self.add_tool("bacgenome.scaffold_abund")  # modified by ghd @ 20190731，binning流程把tool改掉了，旧的tool放到细菌里
        self.coverage2 = self.add_tool("bacgenome.scaf_mean_cov")
        self.correct = self.add_tool("bacgenome.pilon")
        self.correct_reads_is_complete = False

    def check_options(self):
        """
        检查参数
        :return:
        """
        self.prefix = os.path.basename(self.option("fq1").prop["path"]).split(".")[0]
        return True

    def run(self):
        super(MappingModule, self).run()
        self.run_bowtie()

    def run_bowtie(self):
        self.bowtie2.set_options({
            "fastq1": self.option("fq1"),
            "fastq2": self.option("fq2"),
            "ref_fasta": self.option("ref_fa")
        })
        self.bowtie2.on("end", self.run_convert)
        self.bowtie2.run()

    def run_convert(self):
        self.convert.set_options({
            "sam": self.bowtie2.output_dir + "/" + self.prefix + ".pair.sam",
            "analysis": "metagbin"
        })
        if self.option("coverage") and not self.option("cal_cov_auto"):
            self.convert.on("end", self.run_coverage)
            self.convert.on("end", self.run_coverage2)
            self.run_tools.append(self.coverage)
            self.run_tools.append(self.coverage2)
            if self.option("correct_reads"):
                self.correct_reads_is_complete = True
        if self.option("correct_reads"):
            self.convert.on("end", self.run_correct)
            self.run_tools.append(self.correct)
        if len(self.run_tools) == 0:
            self.convert.on("end", self.set_output)
        elif len(self.run_tools) == 1:
            # 此时run_tools里面是self.correct
            if self.option("cal_cov_auto"):
                self.correct.on("end", self.estimate_run_coverage)
            else:
                self.on_rely(self.run_tools, self.set_output)
        else:
            self.on_rely(self.run_tools, self.set_output)
        self.convert.run()

    def run_coverage(self):
        self.coverage.set_options({
            "bam_file": self.convert.output_dir + "/" + self.prefix + "_sorted.bam"
        })
        self.coverage.run()

    def run_coverage2(self):
        self.coverage2.set_options({
            "bam_file": self.convert.output_dir + "/" + self.prefix + "_sorted.bam"
        })
        self.coverage2.run()

    def run_correct(self):
        self.correct.set_options({
            "ref_fa": self.option("ref_fa"),
            "frags_bam": self.convert.output_dir + "/" + self.prefix + "_sorted.bam"
        })
        self.correct.run()

    def estimate_run_coverage(self):
        file = os.path.join(self.correct.work_dir, "pilon.o")
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("Corrected"):
                    # eg.'Corrected 0 snps; 0 ambiguous bases; corrected 57 small insertions totaling 57 bases, 1 small deletions totaling 1 bases'
                    a = line.split(";")  # eg.['Corrected 0 snps', ' 0 ambiguous bases', ' corrected 57 small insertions totaling 57 bases, 1 small deletions totaling 1 bases']
                    snp = int(a[0].split(" ")[1])  # eg.['Corrected', '0', 'snps']
                    ambiguous = int(a[1].split(" ")[1])  # eg. ['', '0', 'ambiguous', 'bases']
                    b = a[2].split(",") # eg. [' corrected 57 small insertions totaling 57 bases', ' 1 small deletions totaling 1 bases']
                    insert = int(b[0].split(" ")[2]) # eg. ['', 'corrected', '57', 'small', 'insertions', 'totaling', '57', 'bases']
                    delete = int(b[1].split(" ")[1]) # eg. ['', '1', 'small', 'deletions', 'totaling', '1', 'bases']
                    # 最初的校正标准是统一校正两次，由杨帆确定，现在去掉了这条标准
                    # if snp > 10 or ambiguous > 10 or insert > 10 or delete > 10:
                    # 上面一行是张超认为校正成功的标准
                    if snp != 0 or ambiguous != 0 or insert != 0 or delete != 0:
                    # 上面一行是李莎莎认为校正成功的标准
                    # 暂时使用李莎莎认为的校正标准，实际上相当于所有的任务统一校正5次
                    # 以后产品线考虑清楚，统一了线下校正标准后再做更新
                        self.set_output()
                        return
        self.correct_reads_is_complete = True
        self.on_rely([self.coverage, self.coverage2], self.set_output)
        self.run_coverage()
        self.run_coverage2()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.option("sam_out", self.bowtie2.output_dir + "/" + self.prefix + ".pair.sam")
        self.option("bam_out", self.convert.output_dir + "/" + self.prefix + "_sorted.bam")
        if self.option("cal_cov_auto") and self.correct_reads_is_complete:
            self.option("cov_out", self.coverage.option("metabat_depth").prop["path"])
            self.option("cov_out2", self.coverage2.option("depth").prop["path"])
        elif self.option("coverage") and not self.option("cal_cov_auto"):
            self.option("cov_out", self.coverage.option("metabat_depth").prop["path"])
            self.option("cov_out2", self.coverage2.option("depth").prop["path"])
        if self.option("correct_reads"):
            self.option("scf_seq", self.correct.option("scf_seq").prop["path"])
        self.logger.info("设置丰度计算结果成功")
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [",", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MappingModule, self).end()
