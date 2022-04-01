# -*- coding: utf-8 -*-
# __author__ = 'xuting'
#last modify qingchen.zhang@20190509
from __future__ import division
import math
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class MgBaseInfoAgent(Agent):
    """
    version 1.0
    author: xuting
    last_modify: 2018.07.12
    用于统计多个fastq文件的碱基质量信息
    """
    def __init__(self, parent):
        super(MgBaseInfoAgent, self).__init__(parent)
        options = [
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq"}]  # 输入fastq文件
        self.add_option(options)
        self.step.add_steps("base_info_stat")
        self.on('start', self.start_base_info)
        self.on('end', self.end_base_info)
        self._memory_increase_step = 20  # 每次重运行增加内存20G by guhaidong @ 20190111

    def start_base_info(self):
        self.step.base_info_stat.start()
        self.step.update()

    def end_base_info(self):
        self.step.base_info_stat.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("fastq_path").is_set:
            raise OptionError("参数fastq不能为空")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        total = os.path.getsize(self.option("fastq_path").prop['path']) / (1024 * 1024 * 1024)
        total = total
        total = math.ceil(total)
        self._memory = '{}G'.format(int(total) + 2)  # 增加2G内存 by qingchen.zhang @ 20190521
        #self._memory = '10G'
        #total = 0
        #for f in self.option("fastq_path").prop["samples"]:
            #total += os.path.getsize(f)
        #total = total / (1024 * 1024 * 1024)
        #total = total * 4
        #total = math.ceil(total)
        #self._memory = '{}G'.format(int(total) + 2)  # 增加2G内存 by ghd @ 20180712

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'base_info', "", "样品碱基信息统计目录"],
            [r".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["\./base_info/.+\.fastxstat\.txt$", "xls", "碱基质量统计文件"]
        ])
        super(MgBaseInfoAgent, self).end()


class MgBaseInfoTool(Tool):
    def __init__(self, config):
        super(MgBaseInfoTool, self).__init__(config)
        self._version = 1.0
        self.fastx_stats_path = "bioinfo/seq/fastx_toolkit_0.0.14/fastx_quality_stats"

    def _run_fastx(self):
        fastq = self.option('fastq_path').prop['path']
        fastq_name = os.path.basename(fastq) + ".fastxstat.txt"
        file_name = os.path.join(self.output_dir, fastq_name)
        cmd = self.fastx_stats_path + " -i " + fastq + " -o " + file_name
        self.logger.info("运行base_info，进行统计")
        command = self.add_command('base_info_cmd', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行base_info完成")
        else:
            self.set_error("ERROR:Sample %s" %(os.path.basename(fastq))) #add by qingchen.zhang @20191012

    def run(self):
        """
        运行
        """
        super(MgBaseInfoTool, self).run()
        self._run_fastx()
        self.end()
