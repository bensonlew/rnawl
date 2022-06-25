# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
import re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class TableSelectAgent(Agent):
    def __init__(self, parent):
        super(TableSelectAgent, self).__init__(parent)
        options = [
            {"name": "origin_table", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "sep", "type": 'string', "default": "tab"},
            {"name": 'select_table', "type": 'outfile', 'format': 'ref_rna_v2.common'},
        ]
        self.add_option(options)
        self.step.add_steps("table_select")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.table_select.start()
        self.step.update()

    def step_finish(self):
        self.step.table_select.finish()
        self.step.update()

    def check_options(self):
        if not self.option("origin_table").is_set:
            raise OptionError("必须设置输入数据表")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(TableSelectAgent, self).end()


class TableSelectTool(Tool):
    def __init__(self, config):
        super(TableSelectTool, self).__init__(config)
        self._version = "v1.0"
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin'))
        self.file = {
            'processed': os.path.join(self.work_dir, 'processed_table.txt'),
        }

    def run(self):
        super(TableSelectTool, self).run()
        self.processing()
        self.set_output()
        self.end()

    def processing(self):
        sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[self.option('sep')]
        df = pd.read_table(self.option('origin_table').prop['path'], header=0, sep=sep)
        df_mean = df.groupby(by=df.columns[0], sort=False, as_index=False).mean()
        df_mean[df.columns[0]] = df_mean[df.columns[0]].apply(lambda x: x.replace(';', '_').replace(',', '_').replace('(', '_').replace(')', '_').replace('\'', '‘'))
        df_mean.to_csv(self.file['processed'], header=True, index=False, sep=sep)

    # def id_processing(self):
    #     sep = {'tab': '\t', 'comma': ',', 'semicolon': ';', 'space': '\s+'}[self.option('sep')]
    #     with open(self.option('origin_table').prop['path'], 'r') as in_f, open(self.file['processed'], 'w') as of:
    #         head = in_f.readline()
    #         of.write(head)
    #         for line in in_f:
    #             items = line.split(sep, 1)
    #             name = items[0].replace(';', '_').replace(',', '_').replace('(', '_').replace(')', '_')
    #             of.write(name + sep + items[1])

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if os.path.exists(self.file['processed']):
            self.logger.info("设置select_table成功")
            self.option("select_table", self.file['processed'])
