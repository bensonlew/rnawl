# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.05.17

import os
import re
import datetime
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess


class ExtractSeqAssembleAgent(Agent):
    """
    细菌基因组用于页面下载序列
    """

    def __init__(self, parent):
        super(ExtractSeqAssembleAgent, self).__init__(parent)
        options = [
            {"name": "seq_path", "type": "string", 'default': ""},
            {"name": "type", "type": "string", 'default': ""},
            {"name": "sample", "type": "string", 'default': ""},
            {"name": "seq_list", "type": "string", 'default': ""},
            {"name": "sample_new", "type": "string", 'default': ""},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("seq_path"):
            raise OptionError("必须设置参数seq_path，提供序列文件路径！", code="31401301")
        if not self.option("sample"):
            raise OptionError("必须设置参数sample！", code="31401302")
        if self.option("type") == "":
            raise OptionError("下载序列类型！", code="31401303")
        if self.option("seq_list") == "":
            raise OptionError("提供需要下载的基因ID！", code="31401304")
        if self.option("sample_new") == "":
            raise OptionError("提供需要下载样品的改名！", code="31401305")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(ExtractSeqAssembleAgent, self).end()


class ExtractSeqAssembleTool(Tool):
    def __init__(self, config):
        super(ExtractSeqAssembleTool, self).__init__(config)

    def run_get_seq(self):
        sample_new = self.option("sample_new")
        sample= self.option("sample")
        seq_list =self.option('seq_list')
        path =self.option('seq_path')
        gene_new = []
        if self.option('type') in ['uncomplete']:
            seqs = seq_list.split(',')
            out = ''
            self.logger.info(seq_list)
            if re.search(r'Scaffold',seq_list):
                out = self.output_dir + '/' + sample_new + '_' + str(len(seqs)) + '_scaffold.fasta'
            elif re.search(r'contig',seq_list):
                out = self.output_dir + '/' + sample_new + '_' + str(len(seqs)) + '_contig.fasta'
            self.logger.info(out)
            lds = []
            for dd in seqs:
                file =path + '/' + dd + '.fasta'
                lds.append(file)
            fil=' '.join(lds)
            cmd = 'cat %s >> %s' %(fil,out)
            try:
                subprocess.check_output(cmd, shell=True)
                self.logger.info('运行cat seq完成')
            except subprocess.CalledProcessError:
                self.set_error('运行cat seq出错', code="31401301")
        elif self.option('type') in ['complete']:
            if os.path.exists(self.output_dir + '/' + sample_new):
                os.removedirs(self.output_dir + '/' + sample_new)
            else:
                os.mkdir(self.output_dir + '/' + sample_new)
            seqs = seq_list.split(',')
            for dd in seqs:
                file = path + '/' + dd + '.fasta'
                new_file =self.output_dir + '/' + sample_new + '/' + sample_new + '_' + dd + '.fasta'
                os.link(file,new_file)

    def run(self):
        super(ExtractSeqAssembleTool, self).run()
        self.run_get_seq()
        self.end()
