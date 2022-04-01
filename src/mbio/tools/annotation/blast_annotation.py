# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2017.04.19
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import re
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.align.blast.blastout_statistics import *


class BlastAnnotationAgent(Agent):
    """
    对blast结果进行evalue, similarity, identity, length, score的筛选，重新对nr、Swissprot进行注释
    """
    def __init__(self, parent):
        super(BlastAnnotationAgent, self).__init__(parent)
        options = [
            {"name": "blastout_table", "type": "infile", "format": "align.blast.blast_table"},
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "score", "type": "float", "default": 0},  # score值
            {"name": "similarity", "type": "float", "default": 0},  # similarity值
            {"name": "identity", "type": "float", "default": 0},  # identity值
        ]
        self.add_option(options)
        self.step.add_steps("blast_anno")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.blast_anno.start()
        self.step.update()

    def step_end(self):
        self.step.blast_anno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("blastout_table").is_set:
            raise OptionError("必须提供BLAST的结果文件")
        if self.option("evalue"):
            if self.option("evalue") > 1e-3:
                raise OptionError("E-value需小于最大值1e-3")
        if self.option("similarity"):
            if self.option("similarity") > 100 or self.option("similarity") < 0:
                raise OptionError("similarity范围为0-100")
        if self.option("identity"):
            if self.option("identity") > 100 or self.option("identity") < 0:
                raise OptionError("identity范围为0-100")
        else:
            pass

    def set_resource(self):
        self._cpu = 10
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [r".*evalue\.xls", "xls", "比对结果E-value分布图"],
            [r".*similar\.xls", "xls", "比对结果相似度分布图"]
        ])
        super(BlastAnnotationAgent, self).end()


class BlastAnnotationTool(Tool):
    def __init__(self, config):
        super(BlastAnnotationTool, self).__init__(config)
        self.python = self.config.SOFTWARE_DIR + '/program/Python/bin/python'
        self.blast_filter = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/blast_filter.py'
        m = re.match(r"(gene_.+)blastout_table.xls", os.path.basename(self.option("blastout_table").prop["path"]))
        n = re.match(r"(.+)blastout_table.xls", os.path.basename(self.option("blastout_table").prop["path"]))
        if m:
            name = m.group(1)
        elif n:
            name = n.group(1)
        else:
            name = ""
        self.blast_path = self.output_dir + "/" + name + "blast.xls"
        self.evalue_path = self.output_dir + "/" + name + "evalue.xls"
        self.similarity_path = self.output_dir + "/" + name + "similar.xls"

    def run(self):
        super(BlastAnnotationTool, self).run()
        self.run_blast_filter()
        self.run_stat()
        self.end()

    def run_stat(self):
        self.logger.info("开始进行统计分析")
        try:
            blastout_statistics(blast_table=self.blast_path, evalue_path=self.evalue_path, similarity_path=self.similarity_path)
            self.logger.info("统计分析完成")
        except Exception as e:
            self.set_error("运行统计出错:{}".format(e))

    def run_blast_filter(self):
        self.logger.info("add_command开始进行blast参数筛选")
        cmd = '{} {} {} {} {} {} {} {}'.format("program/Python/bin/python", self.blast_filter, self.option("blastout_table").prop["path"], self.blast_path, self.option("evalue"), self.option("score"), self.option("similarity"), self.option("identity"))
        command_obj = self.add_command("blast_filter", cmd).run()
        self.wait(command_obj)
        if command_obj.return_code == 0:
            self.logger.info("add_command blast参数筛选完成")
        else:
            self.set_error("add_command blast参数筛选失败")
