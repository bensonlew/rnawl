# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2019.03.21

import os
import re
import datetime
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class ExtractSeqAssAgent(Agent):
    """
    细菌基因组扫描图组装序列下载
    """

    def __init__(self, parent):
        super(ExtractSeqAssAgent, self).__init__(parent)
        options = [
            {"name": "seq_path", "type": "infile","format": "sequence.fasta"},  # 前端传的到样品名的路径
            {"name": "seq_type", "type": "string"},
            {"name": "seq_list", "type": "string"},  # 需要下载的序列list
            {"name": "sample_new", "type": "string"}
            # ‘{"B":"B_new"}第一个是新的样品名
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("seq_path"):
            raise OptionError("必须设置参数seq_path，提供序列文件路径！")
        if not self.option("seq_type"):
            raise OptionError("必须设置参数seq_type！")
        if self.option("seq_list") == "":
            raise OptionError("提供需要下载的基因ID！")
        if self.option("sample_new") == "":
            raise OptionError("提供需要下载的样品sample_new！")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(ExtractSeqAssAgent, self).end()

class ExtractSeqAssTool(Tool):
    def __init__(self, config):
        super(ExtractSeqAssTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.script = self.config.PACKAGE_DIR + '/bacgenome/extract_seq_byass.pl'

    def run_get_seq(self):
        sample_new = self.option("sample_new")
        self.logger.info(sample_new)
        seq_path = self.option("seq_path").prop['path']
        seq_type = self.option("seq_type")
        cmd = "{} {} {} {} {} {} {}".format(self.perl_path, self.script, seq_path,
                                            self.option("seq_list"), sample_new, seq_type, self.output_dir)
        self.logger.info(cmd)
        command = self.add_command("get_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("get_seq succeed")
        else:
            self.set_error("get_seq failed")

    def run(self):
        super(ExtractSeqAssTool, self).run()
        self.run_get_seq()
        self.end()
