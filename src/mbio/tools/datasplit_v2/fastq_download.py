# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modified: 20201105

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class FastqDownloadAgent(Agent):
    """
    下载fastq文件
    """
    def __init__(self, parent=None):
        super(FastqDownloadAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "datasplit.fastq"},
            {"name": "out_fastq", "type": "outfile", "format": "datasplit.fastq"},
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("fastq").is_set:
            raise OptionError("缺少fastq,请设置")

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(FastqDownloadAgent, self).end()


class FastqDownloadTool(Tool):
    def __init__(self, config):
        super(FastqDownloadTool, self).__init__(config)

    def get_fastq(self):
        """
        获取对象存储的fastq文件
        """
        new_fastq = os.path.join(self.output_dir, os.path.basename(self.option("fastq").prop["path"]))
        if os.path.exists(new_fastq):
            os.remove(new_fastq)
        os.link(self.option("fastq").prop["path"], new_fastq)
        self.option("out_fastq", new_fastq)

    def run(self):
        super(FastqDownloadTool, self).run()
        self.get_fastq()
        self.end()
