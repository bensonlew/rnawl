# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
from __future__ import division
import math
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class CazyAlignAgent(Agent):
    """
    version 1.0
    author: zhouxuan
    last_modify: 20170519
    cazy数据库比对
    """
    def __init__(self, parent):
        super(CazyAlignAgent, self).__init__(parent)
        options = [
            {"name": "faa_file", "type": "infile", "format": "sequence.fasta"}]  # 输入fasta文件氨基酸序列
        self.add_option(options)
        self.step.add_steps("cazy_align")
        self.on('start', self.start_cazy_align)
        self.on('end', self.end_cazy_align)

    def start_cazy_align(self):
        self.step.cazy_align.start()
        self.step.update()

    def end_cazy_align(self):
        self.step.cazy_align.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("faa_file").is_set:
            raise OptionError("参数faa_file不能为空")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        self._memory = '20G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [r'base_info', "", "样品碱基信息统计目录"],
        #     [r".", "", "结果输出目录"]
        # ])
        # result_dir.add_regexp_rules([
        #     ["\./base_info/.+\.fastxstat\.txt$", "xls", "碱基质量统计文件"]
        # ])
        super(CazyAlignAgent, self).end()


class CazyAlignTool(Tool):
    def __init__(self, config):
        super(CazyAlignTool, self).__init__(config)
        self._version = 1.0
        self.cazy_align_path = "bioinfo/align/scripts/cazy_align-zx.py"

    def run_cazy_align(self):
        cmd = "{} -i {} -o {}".format(self.cazy_align_path, self.option('faa_file').prop['path'],
                                      self.output_dir + "/gene")
        self.logger.info("start cazy_align")
        command = self.add_command("cazy_align_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cazy_align done")
        else:
            self.set_error("cazy_align error")
            raise Exception("cazy_align error")

    def run(self):
        """
        运行
        """
        super(CazyAlignTool, self).run()
        self.run_cazy_align()
        self.end()
