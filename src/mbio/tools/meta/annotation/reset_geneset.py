# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180708
# last modified : guhaidong

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd
from pymongo import MongoClient
from biocluster.config import Config
import shutil
from biocluster.file import download


class ResetGenesetAgent(Agent):
    """
    宏基因组注释交互分析总览和KEGG pathway基因集创建，合并功能，并生成基因集统计信息
    """

    def __init__(self, parent):
        super(ResetGenesetAgent, self).__init__(parent)
        options = [
            {"name": "gene_length_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "select_genes", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_relative_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "samples", "type": "string", "default": "all"},
            {"name": "geneset_list", "type": "string", "default": "no_merge"},
            {"name": "select_gene_profile", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "overviewfile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "mytype", "type": "int", "choose":[1,2,3,4]},
        ]
        self.add_option(options)
        self.gene_length_table = self.option("gene_length_table")
        self._memory_increase_step = 50  # 每次重运行增加50G内存 add by GHD @ 20180708

    def check_options(self):
        if not self.option("gene_length_table").is_set:
            raise OptionError("必须设置基因长度文件", code="32700901")
        if not self.option("gene_profile").is_set:
            raise OptionError("必须设置基因丰度文件", code="32700902")
        if not self.option("gene_relative_profile").is_set:
            raise OptionError("必须设置基因相对丰度文件", code="32700903")
        '''
        if not self.option("select_genes").is_set:
            raise OptionError("必须设置筛选的基因文件", code="32700904")
        '''
        return True

    def set_resource(self):
        self._cpu = 2
        # self._memory = '10G'  # 改回 by guhaidong @ 20180427
        f_sizes = []
        for f in ["gene_length_table", "gene_profile", "gene_relative_profile",
                  "select_genes", "overviewfile"]:
            if self.option(f).is_set:
                f_size = self.option(f).get_size()
                f_sizes.append(f_size)
        self._memory = str(int(max(f_sizes) / (1024 * 1024 * 1024)) + 5) + 'G'

        # memory = 10 + 10 * self._rerun_time  # 每次重运行增加5G内存 by guhaidong @ 20180417
        # self._memory = "%sG" % memory

    def end(self):
        super(ResetGenesetAgent, self).end()


class ResetGenesetTool(Tool):
    def __init__(self, config):
        super(ResetGenesetTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/annotation/mg_annotation/geneset_reset.py'

    def run(self):
        """
        运行
        :return:
        """
        super(ResetGenesetTool, self).run()
        self.run_geneset()
        self.set_output()
        self.end()

    def run_geneset(self):
        self.logger.info("start geneset select")
        samples = self.option("samples")
        profile = self.option("gene_profile").prop["path"]
        relative_profile = self.option("gene_relative_profile").prop["path"]
        gene_length_table = self.option("gene_length_table").prop["path"]
        overviewfile = self.option("overviewfile").prop["path"]
        mytype = self.option("mytype")
        if self.option("geneset_list") == "no_merge":
            geneset_file_list = self.option("geneset_list")
        else:
            check_list = self.option("geneset_list").split(",")
            geneset_file_list = ",".join(download(i) for i in check_list)
        cmd = self.python_path + ' {} -p {} -rp {} -gl {} -sam {} -m {} -o {} -overview {} -ty {}'.format(self.script,
                                                                                            profile, relative_profile,
                                                                                            gene_length_table, samples,
                                                                                            geneset_file_list,
                                                                                            self.output_dir, overviewfile, mytype)
        if self.option("mytype") != 4:
            gene_file = self.option("select_genes").prop["path"]
            cmd += " -s {} ".format(gene_file)
        command = self.add_command('function_select', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("geneset reset succeed")
        elif command.return_code in [1, -9, -11]:  # add memory limit by guhaidong @ 20180417 add 1 @ 2018.7.25
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("geneset reset failed", code="32700901")

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        try:
            self.option("select_gene_profile", self.output_dir + "/gene_profile/reads_number.xls")
            self.logger.info("设置输出结果文件成功")
        except Exception as e:
            self.set_error("输出结果文件异常——%s", variables=(e), code="32700902")
