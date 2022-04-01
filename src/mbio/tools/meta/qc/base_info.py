# -*- coding: utf-8 -*-
# __author__ = 'xuting'
from __future__ import division
import math
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class BaseInfoAgent(Agent):
    """
    version 1.0
    author: xuting
    last_modify: 2018.07.12
    用于统计多个fastq文件的碱基质量信息
    """
    def __init__(self, parent):
        super(BaseInfoAgent, self).__init__(parent)
        options = [
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"}]  # 输入fastq文件夹
        self.add_option(options)
        self.step.add_steps("base_info_stat")
        self.on('start', self.start_base_info)
        self.on('end', self.end_base_info)

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
            raise OptionError("参数fastq不能为空", code="32706001")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 16
        total = 0
        for f in self.option("fastq_path").prop["samples"]:
            total += os.path.getsize(f)
        total = total / (1024 * 1024 * 1024)
        total = total * 4
        total = math.ceil(total)
        self._memory = '{}G'.format(int(total) + 2)  # 增加2G内存 by ghd @ 20180712

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'base_info', "", "样品碱基信息统计目录"],
            [r".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["\./base_info/.+\.fastxstat\.txt$", "xls", "碱基质量统计文件"]
        ])
        super(BaseInfoAgent, self).end()


class BaseInfoTool(Tool):
    def __init__(self, config):
        super(BaseInfoTool, self).__init__(config)
        self._version = 1.0
        self.fastx_stats_path = "bioinfo/seq/fastx_toolkit_0.0.14/fastx_quality_stats"

    def _run_fastx(self):
        work_path = os.path.join(self.work_dir, "output")
        self.option('fastq_path').get_full_info(work_path)
        base_info_dir = os.path.join(work_path, "base_info")
        if not os.path.exists(base_info_dir):
            os.mkdir(base_info_dir)
        fastq_list = list()
        j = 0
        for fastq in self.option('fastq_path').prop['unzip_fastqs']:
            fastq_list.append(fastq)
            self.logger.info(fastq)
        # 同时运行过多的命令会导致远程的机器失去响应，对同一时间内的运行程序数量做出限制
        while len(fastq_list) > 0:
            cmd_list = list()
            k = 15
            if len(fastq_list) < k:
                k = len(fastq_list)
            for i in range(k):
                j += 1
                self.logger.info("列表长度{}".format(len(fastq_list)))
                fastq = fastq_list.pop()
                file_name = os.path.join(base_info_dir, os.path.basename(fastq) + ".fastxstat.txt")
                cmd = self.fastx_stats_path + " -i " + fastq + " -o " + file_name
                str_ = "fastx_quality_stats_" + str(j) + "_" + re.split("\.", os.path.basename(fastq))[0]
                command = self.add_command(str(str_.lower()), cmd)
                cmd_list.append(command)
            for mycmd in cmd_list:
                self.logger.info("开始运行{}".format(mycmd.name))
                mycmd.run()
            self.wait()
            for mycmd in cmd_list:
                if mycmd.return_code == 0:
                    self.logger.info("运行{}完成".format(mycmd.name))
                else:
                    self.set_error("运行%s出错", variables=(mycmd.name), code="32706001")
                    raise Exception("运行{}出错".format(mycmd.name))

    def run(self):
        """
        运行
        """
        super(BaseInfoTool, self).run()
        self._run_fastx()
        self.end()
