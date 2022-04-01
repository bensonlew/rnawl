# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import numpy as np
from collections import OrderedDict


class StackedColumnAbsoluteAgent(Agent):
    def __init__(self, parent):
        super(StackedColumnAbsoluteAgent, self).__init__(parent)
        options = [
            {"name": "data_file", "type": "infile", "format": "tool_lab.table"},
            {"name": "group_merge", "type": "string", "default": "average"},  # 组内合并，average, or none
            {"name": "group_file", "type": "infile", "format": "tool_lab.group_table"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        data = pd.read_csv(self.option("data_file").prop["path"],'\t')
        data_list = data.columns[1:].tolist()
        list1 = []
        for i in data_list:
            if i not in list1:
                list1 += i
            else:
                raise OptionError('数据文件中样本{}存在重复'.format(i))

        if self.option("group_file").is_set:
            group = pd.read_csv(self.option("group_file").prop["path"], '\t')
            group_list = group.iloc[:, 0].tolist()
            list2 = []
            for i in group_list:
                if i not in list2:
                    list2 += i
                else:
                    raise OptionError('分组文件中样本{}存在重复'.format(i))
            for i in group_list:
                if i not in data_list:
                    raise OptionError('分组文件中的样本{}在数据表中找不到'.format(i))

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(StackedColumnAbsoluteAgent, self).end()


class StackedColumnAbsoluteTool(Tool):
    def __init__(self, config):
        super(StackedColumnAbsoluteTool, self).__init__(config)

    def group_merge(self):  # 组内合并计算
        if self.option("group_merge") == "average":
            with open(self.option("group_file").prop["path"]) as f1:
                dict_group = OrderedDict()
                for line in f1.readlines()[1:]:
                    item = line.strip().split("\t")
                    if item[1] not in dict_group:
                        dict_group[item[1]] = [item[0]]
                    else:
                        dict_group[item[1]].append(item[0])
            f2 = pd.read_table(self.option("data_file").prop["path"], sep='\t')
            name = f2.ix[:, 0]
            data = pd.DataFrame()
            for i in dict_group.keys():
                data[i] = f2.ix[:, dict_group[i]].apply(lambda x: round(x.mean(),2), axis=1) # 根据分组计算平均值
                result_dataframe = pd.concat([name, data], axis=1)
                result_dataframe.to_csv(os.path.join(self.output_dir, "group_result.txt"), index=0, sep='\t',
                                        encoding="utf-8")
        else:  # 组内不合并时，将数据文件直接输出到output_dir中
            output_file = self.option("data_file").prop["path"]   
            output_link = self.output_dir + "/group_result.txt"
            os.link(output_file,output_link)

    def run(self):
        super(StackedColumnAbsoluteTool, self).run()
        self.group_merge()
        self.end()
