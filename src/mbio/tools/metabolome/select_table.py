# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180911
# last modified: guhaidong

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd
import numpy as np
import commands


class SelectTableAgent(Agent):
    """
    根据样本名和gene_list筛选profile并根据样本分组计算组和，
    """

    def __init__(self, parent):
        super(SelectTableAgent, self).__init__(parent)
        options = [
            {"name": "origin_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "select_genes", "type": "infile", "format": "sequence.profile_table"},
            {"name": "samples", "type": "string", "default": "all"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "select_table", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "group_all", "type": "string", "default": "F"},  # group_name 替换为all
            {"name": "group_method", "type": "int", "default": 1},  # 1：求和，2：求均值，3：求中位数
            {"name": "st", "type": "string", "default": "T"},  # 是否计算total
            {"name": "top", "type": "int","default": 0},
            {"name": "select_columns", "type": "string","default": "GeneID"},
            {"name": "merge", "type": "string"},
            {"name": "trans", "type": "bool", "default": False}, #是否转置
            {"name": "scale", "type": "bool", "default": False},  ## 是否标准化,代谢用
            {"name": "select_origin_abu", "type": "outfile", "format": "sequence.profile_table"}, ## 有scale时使用
            {"name": "filter_zero", "type": "string", "default": ""}, #对total为0的丰度进行过滤add by qingchen.zhang@20190315
            {"name": "scale_method","type":"string","default":"UV"}, #zouuguanqing 20190605
            {"name": "log10", "type": "bool", "default": False}  # zhaoyuzhuo 20210119
        ]
        self.add_option(options)
        self._memory_increase_step = 30  # modified by GHD @ 20180628


    def check_options(self):
        if not self.option("origin_table").is_set:
            raise OptionError("必须设置原始文件", code="34002501")
        if not (self.option("select_genes").is_set or self.option("group").is_set or self.option("samples") != "all"):
            raise OptionError("请设置筛选条件基因集或者group分组计算或者挑选样本", code="34002502")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '16G'  #  by guhaidong @ 20180628
        # memory = 8 + 10 * self._rerun_time  # 每次重运行增加5G内存 by guhaidong @ 20180417
        # self._memory = "%sG" % memory

    def end(self):
        super(SelectTableAgent, self).end()


class SelectTableTool(Tool):
    def __init__(self, config):
        super(SelectTableTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/metabolome/scripts/profile_select.py'

    def run(self):
        """
        运行
        :return:
        """
        super(SelectTableTool, self).run()
        self.run_geneset()
        self.set_output()
        self.end()

    def run_geneset(self):
        self.logger.info("start table_select")
        samples = self.option("samples")
        #self.logger.info(self.option("select_genes")._properties.keys())
        #self.logger.info(self.option("select_genes").__dict__)
        outfile = self.output_dir + "/select_table.xls"
        select_columns = self.option("select_columns")
        origin_table = self.option("origin_table").prop["path"]
        df_exp = pd.read_table(origin_table, "\t")
        if self.option("log10"):
            log10_origin_table = self.output_dir + "/exp_log10.xls"
            column_list = df_exp.columns[1::].tolist()
            for i in column_list:
                df_exp[i] = df_exp[i].map(lambda x: 0 if x == 0 else np.log10(x))
            df_exp.to_csv(log10_origin_table, "\t", index=False)
            cmd = self.python_path + ' {} -i {} -o {} -sc {}'.format(self.script, log10_origin_table, outfile, select_columns)
        else:
            cmd = self.python_path + ' {} -i {} -o {} -sc {}'.format(self.script, origin_table, outfile, select_columns)
        if samples != "all":
            cmd += " -sam " + samples
            if "st" in self.get_option_object().keys():
                cmd += " -st " + self.option("st")
        if "top" in self.get_option_object().keys():
            cmd += " -top " + str(self.option("top"))
        if self.option("select_genes").is_set:
            gene_file = self.option("select_genes").prop["path"]
            cmd += " -s " + gene_file
        if self.option("group").is_set:
            self.logger.info("输出group——file！")
            group_file = self.option("group").prop["path"]
            self.logger.info(group_file)
            if self.option("group_all") == "T":
                group_table = pd.read_table(group_file, sep='\t', header=0)
                group_table["group_name"] = group_table["group_name"].replace(".*", "All", regex=True)
                group_file = self.work_dir + "/group_new.xls"
                group_table.to_csv(group_file, sep="\t", index=False)
            cmd += " -g " + group_file + " -gm " + str(self.option("group_method"))
        if self.option("merge"):
            cmd += " -merge " + self.option("merge")
        if self.option("trans"):
            cmd += " -trans True "
        if "scale" in self.get_option_object().keys() and self.option("scale"):
            cmd += " --scale -odir " + self.output_dir
            if 'scale_method' in self.get_option_object().keys():  # zouguanqing 20190605
                cmd += " -scale_method %s" % self.option('scale_method')
        if "filter_zero" in self.get_option_object().keys() and self.option("filter_zero"):
            cmd += " -filter " + self.option("filter_zero")
        self.logger.info(cmd)
        command = self.add_command('table_select', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("table_select succeed")
        elif command.return_code in [-9,1]:  # modified return_code by guhaidong @ 20180628
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("table_select failed", code="34002501")

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        outfile = self.output_dir + "/select_table.xls"
        cmd = '/usr/bin/head -n 10 {} | /usr/bin/wc -l'.format(outfile)
        _, lines = commands.getstatusoutput(cmd)
        if int(lines) > 1:
            self.option("select_table").set_path(outfile)
            self.logger.info("设置输出结果文件成功")
            self.logger.info(self.option("select_table").path)
        else:
            self.set_error("原始表格在该基因集筛选下为空！", code="34002502")
        if "scale" in self.get_option_object().keys() and self.option("scale"):
            before_scale_table = self.output_dir + "/select_before_scale.xls"
            if os.path.exists(before_scale_table):
                self.option("select_origin_abu").set_path(before_scale_table)
