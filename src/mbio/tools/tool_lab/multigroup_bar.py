# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'yifei.fang'@ 20210816
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import numpy as np
import collections
import math

class MultigroupBarAgent(Agent):
    """
    MultigroupBar
    """
    def __init__(self, parent):
        super(MultigroupBarAgent, self).__init__(parent)
        options = [
            {"name": "bar_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "taxonomy_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "combination", "type": "string", "default": "none"},
            {"name": "error_bar", "type": "string", "default": "sd"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("bar_table").is_set:
            raise OptionError("必须输入table文件")
        if not self.option("taxonomy_table").is_set:
            raise OptionError("必须输入taxonomy文件")

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
        super(MultigroupBarAgent, self).end()


class MultigroupBarTool(Tool):
    """
    分类柱形图数据处理
    """
    def __init__(self, config):
        super(MultigroupBarTool, self).__init__(config)

    def run_multigroup_bar(self):
        """
        开始处理数据
        """
        if self.option("combination") != "none":
            self.logger.info("组内合并")
            with open(self.option("bar_table").prop["path"]) as f1, \
                    open(self.option("taxonomy_table").prop["path"]) as f2, \
                       open(os.path.join(self.output_dir, "result_file.txt"), "w") as w1:
                dict_sample = collections.OrderedDict()
                lines = f1.readlines()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    if item[0] not in dict_sample.keys():
                        dict_sample[item[0]] = []
                        a_list = [float(x) for x in item[1: len(item)]]
                        dict_sample[item[0]].append(np.mean(a_list))
                        if self.option("error_bar") == "sd":
                            dict_sample[item[0]].append(np.std(a_list, ddof=1))
                        else:
                            dict_sample[item[0]].append(np.std(a_list, ddof=1)/math.sqrt(len(a_list)))
                    else:
                        print "第一列{}重复，请核实".format(item[0])
                self.logger.info(dict_sample)
                rows = f2.readlines()
                for row in rows[1:]:
                    item1 = row.strip().split("\t")
                    if item1[0] in dict_sample.keys():
                        dict_sample[item1[0]].append(item1[1])
                    else:
                        print "第一列{}样本缺失，请核实".format(item1[0])
                if self.option("error_bar") == "sd":
                    w1.write("name" + "\t" + "mean" + "\t" + "SD" + "\t" + "group" + "\n")
                else:
                    w1.write("name" + "\t" + "mean" + "\t" + "SE" + "\t" + "group" + "\n")
                for i in dict_sample.keys():
                    list1 = dict_sample[i]
                    w1.write(i + "\t" + str(list1[0]) + "\t" + str(list1[1]) + "\t" + str(list1[2]) + "\n")
        else:
            with open(self.option("bar_table").prop["path"]) as f3,  \
                open(self.option("taxonomy_table").prop["path"]) as f4, \
                    open(os.path.join(self.output_dir, "result_file.txt"), "w") as w2:
                dict_sample = collections.OrderedDict()
                lines = f3.readlines()
                rows = f4.readlines()
                for line in lines[1:]:
                    item1 = line.strip().split("\t")
                    if item1[0] not in dict_sample.keys():
                        dict_sample[item1[0]] = []
                        dict_sample[item1[0]].append(item1[1])
                for row in rows[1:]:
                    item2 = row.strip().split("\t")
                    if item2[0] in dict_sample.keys():
                        dict_sample[item2[0]].append(item2[1])
                w2.write("name" + "\t" + "value" + "\t" + "group" + "\n")
                for i in dict_sample.keys():
                    b_list = dict_sample[i]
                    w2.write(i + "\t" + str(b_list[0]) + "\t" + str(b_list[1]) + "\n")

    def run(self):
        """
        运行
        """
        super(MultigroupBarTool, self).run()
        self.run_multigroup_bar()
        self.end()
