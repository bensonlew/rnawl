# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'binbin.zhao'@ 20200617
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import numpy as np
import math


class BarBreakAgent(Agent):
    """
    CorHeatmapAgent:用于生成之间的correlation
    """
    def __init__(self, parent):
        super(BarBreakAgent, self).__init__(parent)
        options = [
            {"name": "bar_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "low_point", "type": "float"},  # 下断点值
            {"name": "high_point", "type": "float"},  # 上断点值
            {"name": "main_id", "type": "string"},
            {"name": "ishape", "type": "string", "default": "sd"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("bar_table"):
            raise OptionError("必须输入snp_table文件")
        if not self.option("low_point"):
            raise OptionError("必须输入下断点值")
        if not self.option("high_point"):
            raise OptionError("必须输入上断点值")
        if not self.option("ishape"):
            raise OptionError("必须输入ishape取值")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        """
        计算结束
        """
        super(BarBreakAgent, self).end()


class BarBreakTool(Tool):
    """
    断轴柱形图绘制
    """
    def __init__(self, config):
        super(BarBreakTool, self).__init__(config)

    def draw_bar_break(self):
        """
        开始计算绘图
        """
        if self.option("group_table").is_set:
            self.logger.info("导表文件设置了")
            with open(self.option("group_table").prop["path"]) as f1, \
                    open(self.option("bar_table").prop["path"]) as f2, \
                     open(os.path.join(self.output_dir, "newFile.txt"), "w") as w:
                dict_group = {}
                dict_sample = {}
                list_group = []
                lines = f1.readlines()
                for line in lines:
                    item = line.strip().split("\t")
                    if item[1] not in dict_group.keys():
                        list_group.append(item[1])
                        dict_group[item[1]] = []
                        dict_group[item[1]].append(item[0])
                    else:
                        dict_group[item[1]].append(item[0])
                self.logger.info(dict_group)
                rows = f2.readlines()
                for row in rows[1:]:
                    item1 = row.strip().split("\t")
                    if item1[0] not in dict_sample.keys():
                        dict_sample[item1[0]] = item1[1]
                    else:
                        print "第一列{}重复，请核实".format(item1[0])
                self.logger.info(dict_sample)
                self.logger.info(list_group)
                w.write("group_name" + "\t" + "mean" + "\t" + "std" + "\n")
                if self.option("ishape") == "sd":
                    for group_name in list_group:
                        list = [float(dict_sample[x]) for x in dict_group[group_name]]
                        w.write(group_name + "\t" + str(round(np.mean(list), 2)) + "\t" + str(round(np.std(list, ddof=1), 2)) + "\n")
                elif self.option("ishape") == "sem":
                    for group_name in list_group:
                        list = [float(dict_sample[x]) for x in dict_group[group_name]]
                        w.write(group_name + "\t" + str(round(np.mean(list), 2)) + "\t" + str(round(np.std(list, ddof=1)/math.sqrt(len(list)), 2)) + "\n")
        else:
            with open(self.option("bar_table").prop["path"]) as f3, \
                    open(os.path.join(self.output_dir, "newFile.txt"), "w") as w2:
                lines = f3.readlines()
                w2.write("sample_name" + "\t" + "value" + "\t" + "std" + "\n")
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    w2.write(item[0] + "\t" + item[1] + "\t" + "NA" + "\n")

    def set_db(self):
        self.logger.info("开始导表")
        file_path = os.path.join(self.output_dir, "newFile.txt")
        api_bar_break = self.api.api('tool_lab.bar_break')
        api_bar_break.add_bar_break(self.option('main_id'), file_path, self.option('low_point'), self.option('high_point'))

    def run(self):
        """
        运行
        """
        super(BarBreakTool, self).run()
        self.draw_bar_break()
        self.set_db()
        self.end()
