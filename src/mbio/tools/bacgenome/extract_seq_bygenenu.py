# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# version 1.0
# last_modify: 2018.03.26

import os
import re
import datetime
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class ExtractSeqBygenenuAgent(Agent):
    """
    细菌基因组用于页面下载序列
    """

    def __init__(self, parent):
        super(ExtractSeqBygenenuAgent, self).__init__(parent)
        options = [
            {"name": "seq_prefix_path", "type": "string", 'default': ""},  # 前端传的到样品名的路径
            {"name": "type", "type": "string", 'default': ""},  # gene,trna,rrna
            {"name": "sample", "type": "string", 'default': ""},
            {"name": "gene_list", "type": "string", 'default': ""},  # 需要下载的基因序列list
            {"name": "sample_gene_new", "type": "string", 'default': "\"\""}
            # ‘{"B":"B_new","gene": "gene_new"}’第一个是新的样品名，后面是基因前缀改名
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("seq_prefix_path"):
            raise OptionError("必须设置参数seq_prefix_path，提供序列文件路径！", code="31401401")
        if not self.option("sample"):
            raise OptionError("必须设置参数sample！", code="31401402")
        if self.option("type") == "":
            raise OptionError("下载序列类型！", code="31401403")
        if self.option("gene_list") == "":
            raise OptionError("提供需要下载的基因ID！", code="31401404")
        if self.option("sample_gene_new") == "":
            raise OptionError("提供需要下载的样品和基因的改名对应关系！", code="31401405")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(ExtractSeqBygenenuAgent, self).end()


class ExtractSeqBygenenuTool(Tool):
    def __init__(self, config):
        super(ExtractSeqBygenenuTool, self).__init__(config)
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.script = self.config.PACKAGE_DIR + '/bacgenome/extract_seq_bygenenu.pl'

    def run_get_seq(self):
        sample_gene_new = eval(self.option("sample_gene_new"))
        self.logger.info(sample_gene_new)
        sample_new = list(set(sample_gene_new.keys()))[0]
        self.logger.info(sample_new)
        gene_new = []
        for ge in sample_gene_new.keys():  # 修改基因前缀的功能
            new = ge + ":" + sample_gene_new[ge]
            gene_new.append(new)
        gene_new = ','.join(gene_new)
        seq_path =self.option("seq_prefix_path")
        cmd = "{} {} {} {} {} {} {} {}".format(self.perl_path, self.script, seq_path, self.option("type"),
                                               self.option("gene_list"), sample_new, gene_new, self.output_dir)
        self.logger.info(cmd)
        command = self.add_command("get_gene_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("get_gene_seq succeed")
        else:
            self.set_error("get_gene_seq failed", code="31401401")

    def run(self):
        super(ExtractSeqBygenenuTool, self).run()
        self.run_get_seq()
        self.end()
