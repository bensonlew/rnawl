# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os,json
import linecache
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SeqStatAgent(Agent):
    """
    传入一个group表，以及是否进行样本合并的参数生成一张丰度表并对并依照group表OTU表进行筛选合并
    """
    def __init__(self, parent):
        super(SeqStatAgent, self).__init__(parent)
        options = [
            {"name": "total_genome", "type": "infile", "format": "sequence.profile_table"},  # 对应的参数列表
            {"name": "sequence_dir", "type": "string"},  # 序列seq的dir
            {"name": "raw_dir", "type": "string"},  # 序列的total的dir
        ]
        self.add_option(options)
        self.step.add_steps("seq_stat")
        self.on('start', self.start_seq_stat)
        self.on('end', self.end_seq_stat)

    def start_seq_stat(self):
        self.step.seq_stat.start()
        self.step.update()

    def end_seq_stat(self):
        self.step.seq_stat.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option("total_genome").is_set:
            raise OptionError("必须提供输入的参数列表")

    def end(self):
        super(SeqStatAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "10G"


class SeqStatTool(Tool):
    def __init__(self, config):
        super(SeqStatTool, self).__init__(config)
        self.logger.info("开始进行seq_stat的序列统计")
        self.stat_seq = os.path.join(self.config.PACKAGE_DIR, "bac_comp_genome/seq_stat.py")
        self.python = "/program/Python/bin/python"

    def seq_stat(self):
        """
        统计基因组水平的基础信息
        :return:
        """
        self.logger.info("开始用package进行计算统计")
        out_file = os.path.join(self.work_dir, "upload_file.xls")
        cmd = "{} {} -i {} -seq {} -total {} -o {}".format(self.python, self.stat_seq, self.option("total_genome").prop["path"], self.option("sequence_dir"), self.option("raw_dir"), out_file)
        self.logger.info(cmd)
        command_name = "seq_stat"
        command = self.add_command(command_name, cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("统计序列成功")
            self.set_output()
        else:
            self.set_error("统计序列失败!")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在生成结果文件目录")
        outfile = os.path.join(self.output_dir, "upload_file.xls")
        if os.path.exists(outfile):
            os.remove(outfile)
        os.link(os.path.join(self.work_dir, "upload_file.xls"), outfile)
        self.logger.info("完成结果文件的设置")

    def run(self):
        super(SeqStatTool, self).run()
        self.logger.info("开始运行seq_stat的tool！")
        self.seq_stat()
        self.set_output()
        self.end()
