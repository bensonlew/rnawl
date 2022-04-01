# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'yifei.fang'@ 20210818
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import numpy as np
import pandas as pd
import collections

class AbundanceBubbleAgent(Agent):
    """
    AbundanceBubbleAgent
    """
    def __init__(self, parent):
        super(AbundanceBubbleAgent, self).__init__(parent)
        options = [
            {"name": "data_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "taxonomy_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "sample_label", "type": "string", "default": "col"},
            {"name": "combination", "type": "string", "default": "none"},
            {"name": "taxonomy", "type": "string", "default": "none"},
            {"name": "normalize", "type": "string", "default": "none"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("data_table").is_set:
            raise OptionError("必须输入table文件")


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
        super(AbundanceBubbleAgent, self).end()


class AbundanceBubbleTool(Tool):
    """
    气泡丰度计算
    """
    def __init__(self, config):
        super(AbundanceBubbleTool, self).__init__(config)

    def formattable(self, tablepath):
        """
        转置表格
        """
        this_table = tablepath
        if self.option('sample_label') != 'col':
            newtable = this_table + '.T'
            self.t_table(this_table, newtable)
            return newtable
        else:
            return this_table

    def t_table(self, table_file, new_table):
        """
        转置表格实现
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)

    def calculate_bubble_abundance(self):
        """
        开始计算气泡丰度
        """
        if self.option("group_table").is_set:
            self.logger.info("分组文件设置了")
            with open(self.option("group_table").prop["path"]) as f1:
                dict_group = collections.OrderedDict()
                list_group = []
                lines = f1.readlines()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    if item[1] not in dict_group.keys():
                        list_group.append(item[1])
                        dict_group[item[1]] = []
                        dict_group[item[1]].append(item[0])
                    else:
                        dict_group[item[1]].append(item[0])

            f2 = pd.read_table(self.formattable(self.option("data_table").prop["path"]), sep="\t")
            ID = f2.ix[:, 0]
            data = pd.DataFrame()
            for i in dict_group.keys():
                data[i] = f2.ix[:, dict_group[i]].apply(lambda x: x.mean(), axis=1)
            group_result = pd.concat([ID, data], axis=1)
            group_result.to_csv(os.path.join(self.output_dir, "group_result.txt"), sep='\t', index=False)

            if self.option("normalize") == "z_score":
                normalized_result = pd.concat([ID, (data - data.mean()) / data.std()], axis=1)
            elif self.option("normalize") == "log10":
                normalized_result = pd.concat([ID, np.log10(data)], axis=1)
            elif self.option("normalize") == "relative":
                normalized_result = pd.concat([ID, data / data.sum(numeric_only=True).sum()], axis=1)
            else:
                normalized_result = pd.concat([ID, data], axis=1)
            normalized_result.to_csv(os.path.join(self.output_dir, "normalized_result.txt"), sep="\t", index=0)

        else:
            f3 = pd.read_csv(self.formattable(self.option("data_table").prop["path"]), sep="\t")
            ID = f3.iloc[:, 0]
            data = f3.iloc[:, 1:]
            if self.option("normalize") == "z_score":
                normalized_result = pd.concat([ID, (data - data.mean()) / data.std()], axis=1)
            elif self.option("normalize") == "log10":
                normalized_result = pd.concat([ID, np.log10(data)], axis=1)
            elif self.option("normalize") == "relative":
                normalized_result = pd.concat([ID, data / data.sum(numeric_only=True).sum()], axis=1)
            else:
                normalized_result = pd.concat([ID, data], axis=1)
            normalized_result.to_csv(os.path.join(self.output_dir, "normalized_result.txt"), sep="\t", index=0)

    def run(self):
        """
        运行
        """
        super(AbundanceBubbleTool, self).run()
        self.calculate_bubble_abundance()
        self.end()
