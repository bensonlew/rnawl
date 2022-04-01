# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180706
# last modified : guhaidong

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd
from pymongo import MongoClient
from biocluster.config import Config
import shutil
from mbio.packages.ref_rna.trans_step import step_count
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
import json


class OverviewSelectAgent(Agent):
    """
    宏基因组注释交互分析总览基因集创建部分基因选择
    """

    def __init__(self, parent):
        super(OverviewSelectAgent, self).__init__(parent)
        options = [
            # {"name": "anno_overview_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "select", "type": "string"},
            {"name": "select_type", "type": "string", "default": "2"}, # 1: 或；2:与
            {"name": "type", "type": "int", "default": 1,"choose":[1,2,3,4]},  # 1：注释总览创建基因集，2：kegg部分创建基因集
            {"name": "database_list", "type": "string", "default": "cog,genus,phylum,kegg"},  # type为1时使用
            {"name": "overview_file", "type": "infile", "format": "sequence.profile_table"},  # 是否从线下anno_overview筛选
            {"name": "gene_kegg_anno", "type": "infile", "format": "sequence.profile_table"},  #type为2时使用
            {"name": "gene_list", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "anno_overview", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self._project_type = "metagenomic"
        self._memory_increase_step = 50  # 每次重运行增加50G内存 add by GHD @ 20180706

    def check_options(self):
        if not self.option("select"):
            raise OptionError("必须设置筛选条件", code="32700701")
        '''
        if not (self.option("type") == 1 or self.option("type") == 2 or self.option("type") == 3 or self.option("type") == 4):
            raise OptionError("type类型必须为1或2", code="32700702")
        '''
        if not self.option("gene_kegg_anno").is_set and self.option("type") == 2:
            raise OptionError("type类型为2时必须输入gene_kegg_file文件", code="32700703")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '50G'  # 改回 by guhaidong @ 20180427
        # memory = 10 + 10 * self._rerun_time  # 每次重运行增加5G内存 by guhaidong @ 20180417
        # self._memory = "%sG" % memory

    def end(self):
        super(OverviewSelectAgent, self).end()


class OverviewSelectTool(Tool):
    def __init__(self, config):
        super(OverviewSelectTool, self).__init__(config)
        self.python_path = "/program/Python/bin/python"
        self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/overview.py'

    def run(self):
        """
        运行
        :return:
        """
        super(OverviewSelectTool, self).run()
        self.run_overview_select()
        self.set_output()
        self.end()

    def run_overview_select(self):
        self.logger.info("start overview select")
        select = json.dumps(self.option("select"))
        select = select.replace("$","\$")
        self.logger.info(select)
        cmd = self.python_path + ' {} -s {} -type {} -o {}'.format(self.script, select, self.option("type"),
                                                                   self.output_dir)
        if self.option("type") == 1:
            cmd += " -database " + self.option("database_list")
        elif self.option("type") == 2:
            gene_file = self.option("gene_kegg_anno").prop["path"]
            cmd += " -kegg " + gene_file
        if self.option("overview_file").is_set:
            cmd += " -of " + self.option("overview_file").prop["path"]
        if self.option("task_id"):
            cmd += " -id " + self.option("task_id")
        if self.option("select_type"):
            cmd += " -st " + self.option("select_type")
        if self.config.DBVersion:
            cmd += ' -mongo ' + str(self.config.DBVersion)
        command = self.add_command('overview_select', cmd, ignore_error=True).run()
        self.wait(command)
        err_msg = "根据条件筛选之后的基因集为空，请重新填写筛选的条件！"
        if command.return_code == 0:
            self.logger.info("overview_select succeed")
        elif command.return_code in [-9]:  # add memory limit by shaohua.yuan @ 20180428 add -9 code by haidong.gu @ 20180706
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error(err_msg)

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        try:
            self.option("gene_list", self.output_dir + "/gene_list.xls")
            self.logger.info("设置gene_list结果文件成功")
            if self.option("type") == 1:
                self.option("anno_overview", self.output_dir + "/anno_overview.xls")
            elif self.option("type") == 2:
                self.option("anno_overview", self.output_dir + "/anno_kegg.xls")
            self.logger.info("设置anno_overvier结果文件成功")
        except Exception as e:
            self.set_error("输出结果文件异常——%s", variables=(e), code="32700702")
