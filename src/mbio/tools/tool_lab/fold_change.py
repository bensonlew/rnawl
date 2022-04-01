# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'yifei.fang'@ 20210817
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import numpy as np
import pandas as pd



class FoldChangeAgent(Agent):
    """
    FoldChangeAgent
    """
    def __init__(self, parent):
        super(FoldChangeAgent, self).__init__(parent)
        options = [
            {"name": "bar_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "sample_label", "type": "string", "default": "col"},
            {"name": "control_group", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("bar_table").is_set:
            raise OptionError("必须输入table文件")
        if not self.option("group_table").is_set:
            raise OptionError("必须输入group文件")
        if not self.option("control_group"):
            raise OptionError("必须输入对照组")


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
        super(FoldChangeAgent, self).end()


class FoldChangeTool(Tool):
    """
    差异倍数计算
    """
    def __init__(self, config):
        super(FoldChangeTool, self).__init__(config)

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

    def calculate_fold_change(self):
        """
        开始计算差异倍数
        """
        with open(self.option("group_table").prop["path"], "r") as f1:
            dict_group = {}
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

        f2 = pd.read_table(self.formattable(self.option("bar_table").prop["path"]), sep="\t")
        ID = f2.ix[:, 0]
        data = pd.DataFrame()
        for i in dict_group.keys():
            data[i] = f2.ix[:, dict_group[i]].apply(lambda x: x.mean(), axis=1)
            result_dataframe = pd.concat([ID, data], axis=1)
        if self.option("control_group") == list_group[0]:
            result_dataframe["FC"] = result_dataframe[list_group[0]] / result_dataframe[list_group[1]]
        else:
            result_dataframe["FC"] = result_dataframe[list_group[1]] / result_dataframe[list_group[0]]

        result_dataframe["log2(FC)"] = result_dataframe["FC"].apply(np.log2)
        result = result_dataframe.drop(result_dataframe.iloc[:, [1, 2]], axis=1)
        result.to_csv(os.path.join(self.output_dir, "result.txt"), sep='\t', index=0)

    def run(self):
        """
        运行
        """
        super(FoldChangeTool, self).run()
        self.calculate_fold_change()
        self.end()
