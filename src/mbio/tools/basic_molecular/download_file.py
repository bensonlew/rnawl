# -*- coding: utf-8 -*-
# __author__ = "薛钦文"
# last_modified: 202112077

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class DownloadFileAgent(Agent):
    """
    下载fastq文件
    """
    def __init__(self, parent=None):
        super(DownloadFileAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "out_path", "type": "string"},
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("input_file").is_set:
            raise OptionError("缺少input_file,请设置")
        if not self.option("out_path"):
            raise OptionError("输出文件的路径需要指定")

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(DownloadFileAgent, self).end()


class DownloadFileTool(Tool):
    def __init__(self, config):
        super(DownloadFileTool, self).__init__(config)

    def get_file(self):
        """
        获取对象存储文件
        """
        new_file = os.path.join(self.option("out_path"), os.path.basename(self.option("input_file").prop["path"]))
        if os.path.exists(new_file):
            os.remove(new_file)
        os.link(self.option("input_file").prop["path"], new_file)
        

    def run(self):
        super(DownloadFileTool, self).run()
        self.get_file()
        self.end()
