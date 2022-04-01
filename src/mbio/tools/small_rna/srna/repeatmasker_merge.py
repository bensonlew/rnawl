# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
# version 1.0

import os
import re
import glob
import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest


class RepeatmaskerMergeAgent(Agent):
    """
    RepeatModeler 构建denovo数据库
    RepeatMasker 查找散在重复序列(interspersed repeats)
    TRF 寻找串联重复序列 (tandem repeat)
    """

    def __init__(self, parent):
        super(RepeatmaskerMergeAgent, self).__init__(parent)
        options = [
            {"name": "repeatmasker_merge", "type": "string", "default": ""},  # 合并的repeatmasker结果
            {"name": "input_genome", "type": "infile", "format": "small_rna.fasta"},  # 组装拼接好的scaffold文件
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("repeatmasker_merge") == None:
            raise OptionError("必须设置参数repeatmasker_merge")

    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        super(RepeatmaskerMergeAgent, self).end()


class RepeatmaskerMergeTool(Tool):
    def __init__(self, config):
        super(RepeatmaskerMergeTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/small_rna/"
        self.repeatmasker_path = "/bioinfo/Genomic/Sofware/RepeatMasker/"
        self.cufflinks_path = '/bioinfo/rna/cufflinks-2.2.1/'
        self.bedtools = "/bioinfo/rna/bedtools2-master/bin/"

    def run_repeat_to_gff(self):
        file = self.option("repeatmasker_merge")
        merge_file = os.path.join(self.work_dir, os.path.basename(file))
        if os.path.exists(merge_file):
            os.remove(file)
        os.link(file, merge_file)
        cmd = '{} {}repeat_to_gff.pl {}'.format(self.perl_path, self.perl_script, merge_file)
        self.logger.info(cmd)
        command = self.add_command("run_repeat_to_gff", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("生成gff文件运行完成")
            self.run_extract_seq()
        else:
            self.set_error("生成gff文件运行出错!")

    def run_extract_seq(self):
        """
        根据repeatmasker的输出文件提取fasta序列
        """
        gff = os.path.join(self.work_dir, os.path.basename(self.option("repeatmasker_merge")) + ".gff")
        index = os.path.join(self.work_dir, os.path.basename(self.option("repeatmasker_merge")) + ".bed")
        with open (gff, "r") as f, open (index, "w") as w:
            f.readline()
            for line in f:
                items = line.strip().split("\t")
                name = items[8].split(";")[0].split("=")[1]
                item = [items[0], str(int(items[3])-1), items[4], name, items[5], items[6]]
                w.write("\t".join(item) + "\n")
        cmd = self.bedtools + "bedtools getfasta -fi %s -bed %s -fo %s -name" % (self.option("input_genome").prop["path"], index, os.path.join(self.work_dir, os.path.basename(self.option("repeatmasker_merge")) + ".fa"))
        command = self.add_command("extract_fasta", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fasta序列提取完成")
            self.set_output()
        else:
            self.set_error("fasta序列提取出错!")

    def set_output(self):
        gff = os.path.join(self.work_dir, os.path.basename(self.option("repeatmasker_merge")) + ".gff")
        repeat = os.path.join(self.work_dir, os.path.basename(self.option("repeatmasker_merge")) + ".fa")
        gff_path = self.output_dir + "/repeatmasker_merge.SSR.gff"
        repeat_path = self.output_dir + "/repeatmasker_merge.repeat.fa"
        if os.path.exists(gff_path):
            os.remove(gff_path)
        os.link(gff, gff_path)
        if os.path.exists(repeat_path):
            os.remove(repeat_path)
        os.link(repeat, repeat_path)


    def run(self):
        super(RepeatmaskerMergeTool, self).run()
        self.run_repeat_to_gff()
        self.end()
