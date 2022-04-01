# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
# import subprocess
import glob


class SamtoolsAgent(Agent):
    """
    samtools:sam格式文件处理工具集
    bcftools：bcf
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11
    """

    Tools = ["dict", "faidx", "index", "mpileup", "sort", "view", "call"]

    def __init__(self, parent):
        super(SamtoolsAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列
            {"name": "sam", "type": "infile", "format": "align.bwa.sam"},    # sam格式文件
            {"name": "method", "type": "string", "default": ""},     # samtool工具
            {"name": "in_bam", "type": "infile", "format": "align.bwa.bam"},     # bam格式输入文件(sam的二进制文件)
            {"name": "mpileup_out", "type": "string", "default": "pileup"},  # mpileup 输出格式
            # {"name": "vcf", "type": "outfile", "format": "vcf"},     # Variant Call Format
            {"name": "pileup", "type": "outfile", "format": "denovo_rna.gene_structure.pileup"},  # pileup格式文件
            # {"name": "bcf", "type": "outfile", "format": "bcf"},     # Variant Call Format
            {"name": "out_bam", "type": "outfile", "format": "align.bwa.bam"}  # bam格式输入文件
        ]
        self.add_option(options)
        self.step.add_steps('samtools')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.samtools.start()
        self.step.update()

    def step_end(self):
        self.step.samtools.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if self.option("method") not in self.Tools:
            raise OptionError("请选择正确的工具类型")
        if self.option("method") in ["dict", "faidx", "mpileup"] and not self.option("ref_fasta").is_set:
            raise OptionError("缺少参考序列")
        # if self.option("method") in ["index", "sort"] and not self.option("in_bam").is_set:
        #     raise OptionError("请传入bam格式文件")
        if self.option("method") in ["view", "index", "sort", "mpileup"]:
            if not self.option("sam").is_set and not self.option("in_bam").is_set:
                raise OptionError("请传入sam或bam文件")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '2G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            # ["./estimators.xls", "xls", "alpha多样性指数表"]
        ])
        # print self.get_upload_files()
        super(SamtoolsAgent, self).end()


class SamtoolsTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(SamtoolsTool, self).__init__(config)
        self.samtools_path = "bioinfo/align/samtools-1.3.1/"
        # self.bcftools_path = "rna/bcftools-1.3.1/"
        if self.option("ref_fasta").is_set:
            self.fasta_name = self.option("ref_fasta").prop["path"].split("/")[-1]
        if self.option("sam").is_set:
            self.sam_name = self.option("sam").prop["path"].split("/")[-1]
        if self.option("in_bam").is_set:
            self.bam_name = self.option("in_bam").prop["path"].split("/")[-1]

    def dict(self):
        cmd = "{}samtools dict {} -o {}".format(self.samtools_path, self.option("ref_fasta").prop["path"], self.fasta_name + ".dict")
        # print cmd
        self.logger.info("开始运行samtools dict命令")
        command = self.add_command("dict", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("构建dict完成！")
        else:
            self.set_error("构建dict过程出错")

    def faidx(self):
        cmd = "{}samtools faidx {}".format(self.samtools_path, self.option("ref_fasta").prop["path"])
        # print cmd
        self.logger.info("开始运行samtools faidx命令")
        command = self.add_command("faidx", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("faidx运行完成！")
        else:
            self.set_error("faidx运行过程出错")

    def index(self, bam):
        cmd = "{}samtools index {}".format(self.samtools_path, bam)
        # print cmd
        self.logger.info("开始运行samtools index命令")
        command = self.add_command("index", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("index运行完成！")
        else:
            self.set_error("index运行过程出错")

    def view(self):
        # self.logger.info("run3")
        if not self.option("ref_fasta").is_set:
            ref_command = ""
            ref_t = ""
        else:
            ref_command = self.option("ref_fasta").prop["path"] + ".fai"
            ref_t = "t"
        cmd = "{}samtools view -b{} {} -o {} {}".format(self.samtools_path, ref_t,  ref_command, self.sam_name + ".bam",
                                                        self.option("sam").prop["path"])
        self.logger.info(cmd)
        # print cmd
        self.logger.info("开始运行samtools view命令")
        command = self.add_command("view", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("view运行完成！")
        else:
            self.set_error("view运行过程出错")

    def sort(self, bam, sorted_bam):
        cmd = "{}samtools sort -o {} {}".format(self.samtools_path, sorted_bam, bam)
        # print cmd
        self.logger.info("开始运行samtools sort命令")
        command = self.add_command("sort", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行sort完成！")
        else:
            self.set_error("运行sort过程出错")

    def mpileup(self, mpileup_out, sort_bam):
        out_option = {"pileup": "", "bcf": "-gu", "vcf": "-vu"}
        cmd = "{}samtools mpileup {} -f {} -o {} {}".format(self.samtools_path, out_option[self.option("mpileup_out")], self.option("ref_fasta").prop["path"], mpileup_out,  sort_bam)
        # print cmd
        self.logger.info(cmd)
        self.logger.info("开始运行samtools mpileup命令")
        command = self.add_command("mpileup", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info(" mpileup命令完成！")
        else:
            self.set_error(" mpileup命令过程出错")

    def set_ouput(self):
        self.logger.info("set out put")
        postfix = ""
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        if self.option("method") in ["index", "sort"]:
            postfix = "bam"
        elif self.option("method") in ["mpileup"]:
            postfix = "pileup"
        file_path = glob.glob(r"*.{}".format(postfix))
        for f in file_path:
            output_dir = os.path.join(self.output_dir, f)
            os.link(os.path.join(self.work_dir, f), output_dir)
            # os.remove(os.path.join(self.work_dir, f))
            if postfix == "out_bam":
                self.option("out_bam").set_path(output_dir)
            else:
                self.option("pileup").set_path(output_dir)
        self.logger.info("output done")
        # self.end()

    def run(self):
        """
        运行
        """
        super(SamtoolsTool, self).run()
        # self.logger.info("run1")
        if self.option("ref_fasta").is_set and not self.option("method") in ["dict", "faidx"]:
            self.faidx()
        self.logger.info(self.option("method"))
        if self.option("method") == "dict":
            self.dict()
        elif self.option("method") == "faidx":
            self.faidx()
        elif self.option("method") == "index" and self.option("sam").is_set:
            self.view()
            self.sort(self.sam_name + ".bam", self.sam_name + ".sorted" + ".bam")
            self.index(self.sam_name + ".sorted" + ".bam")
        elif self.option("method") == "index" and self.option("in_bam").is_set:
            self.index(self.option("in_bam").prop["path"])
        elif self.option("method") == "view":
            self.view()
        elif self.option("method") == "sort" and self.option("sam").is_set:
            self.view()
            self.sort(self.sam_name + ".bam", self.sam_name + ".bam.sorted.bam")
        elif self.option("method") == "sort" and self.option("in_bam").is_set:
            self.sort(self.option("in_bam").prop["path"], self.bam_name + "sorted.bam")
        elif self.option("method") == "mpileup" and self.option("sam").is_set:
            self.view()
            self.sort(self.sam_name + ".bam", self.sam_name + ".sorted.bam")
            self.index(self.sam_name + ".sorted.bam")
            self.mpileup(self.sam_name + ".sorted.bam.plileup", self.sam_name + ".sorted.bam")
        elif self.option("method") == "mpileup" and self.option("in_bam").is_set:
            # self.sort(self.sam_name + ".bam", self.sam_name + ".sorted.bam")
            self.index(self.option("in_bam").prop["path"])
            self.mpileup(self.bam_name + ".pileup", self.option("in_bam").prop["path"])
        else:
            self.logger.info("False")
        self.set_ouput()
        self.end()
