# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# last_modify: 20181218

import os
import re
import math
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SstacksAgent(Agent):
    """
    fastq均一化
    """
    def __init__(self, parent=None):
        super(SstacksAgent, self).__init__(parent)
        options = [
            {"name": "sample_name", "type": "string"},
            {"name": "catlog_dir", "type": "infile", "format": "noref_wgs.catlog_dir"},  # 需要修改infile文件
            {"name": "sample_loci_dir", "type": "infile", "format": "noref_wgs.sample_loci_dir"},  # 需要修改infile文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sample_name"):
            raise OptionError("请输入样品名称", code="35500305")
        if not self.option("catlog_dir"):
            raise OptionError("请输入catlog_dir文件夹", code="35500306")

    def set_resource(self):
        self._cpu = 17   # 因为下面线程为16
        self._memory = "100G"

    def end(self):
        super(SstacksAgent, self).end()


class SstacksTool(Tool):
    def __init__(self, config):
        super(SstacksTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64")
        self.sstacks = "/bioinfo/noRefWGS/stacks-2.2/bin/sstacks"

    def run_sstacks(self):

        cmd = "{} -s {} -c {} -p 16 -o {}"\
            .format(self.sstacks, self.option("sample_loci_dir").prop["path"] + "/" + self.option("sample_name"),
                    self.option("catlog_dir").prop["path"], self.output_dir)
        command = self.add_command("sstacks", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("sstacks运行完成")
        else:
            self.set_error("sstacks运行失败", code="35500303")

    def run(self):
        super(SstacksTool, self).run()
        self.run_sstacks()
        self.end()
