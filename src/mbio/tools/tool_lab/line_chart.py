# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

import os
import pandas as pd
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class LineChartAgent(Agent):
    def __init__(self, parent):
        super(LineChartAgent, self).__init__(parent)
        options = [
            {"name": "data_file", "type": "infile", "format": "tool_lab.table"},  # 输入数据表
            {"name": "ranks", "type": "string", "default": "row"},  # 选择行列标签, row or column
            {"name": "line_type", "type": "string", "default": ''},  # 选择单折线或多折线，single or multi
            {"name": "group_file", "type": "infile", "format": "tool_lab.group_table"},  # 输入分组文件
            {"name": "group_repetition", "type": "string", "default": "true"}  # 组内重复, true or none
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(LineChartAgent, self).end()


class LineChartTool(Tool):
    def __init__(self, config):
        super(LineChartTool, self).__init__(config)

    def t_table(self, table_file, new_table):
        """
        转置表格实现
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)

    def file_check(self):
        if self.option('ranks') == "row":
            data_table = self.option('data_file').prop['path']
        else:
            data_table = self.output_dir + "/data_file_t.txt"
            self.t_table(self.option('data_file').prop['path'], data_table)
        df_data = pd.read_csv(data_table, sep='\t')
        df_group = pd.read_csv(self.option("group_file").prop['path'], sep='\t')
        sample_list1 = df_data.iloc[:, 0].tolist()
        print(sample_list1)
        sample_list2 = df_group.iloc[:, 0].tolist()
        list1 = []
        list2 = []
        for i in sample_list1:
            if i not in list1:
                list1.append(i)
            else:
                raise OptionError("数据表中样本{}存在重复，请检查后重新上传".format(i))
        for i in sample_list2:
            if i not in list2:
                list2.append(i)
            else:
                raise OptionError("分组文件中样本{}存在重复，请检查后重新上传".format(i))
        for a in sample_list2:  # 检查分组样本和数据表中样本是否一致
            if a not in sample_list1:
                raise OptionError("分组文件中的样本{}在数据表中找不到".format(a))
        if self.option("line_type") == 'single':  # 单折线图
            df_count = df_group.iloc[:, 0].groupby(df_group.iloc[:, 1]).count().reset_index(drop=False)
            for i in df_count.index:
                if self.option("group_repetition") == 'true':  # 如果选择组内重复
                    if df_count.iloc[i, 1] < 3:
                        raise OptionError("分组{}中至少应包含三个样本".format(df_count.iloc[i, 0]))
                else:  # 如果选择组内不重复
                    if df_count.iloc[i, 1] != 1:
                        raise OptionError("分组{}中存在重复".format(df_count.iloc[i, 0]))
        else:  # 多折线图
            if len(df_group.columns) != 3:
                raise OptionError('分组文件必须为三列')
            df_count = df_group.iloc[:, 0].groupby([df_group.iloc[:, 1], df_group.iloc[:, 2]]).count().reset_index(
                drop=False)  # 计算分组中样本数
            list_group = df_count.iloc[:, 0].tolist()
            list_subgroup = df_count.iloc[:, 1].tolist()
            if (len(set(list_group)) or len(set(list_subgroup))) < 2:  # 大组和亚组的分组数不小于2
                raise OptionError('分组文件不满足参数要求')
            for i in df_count.index:
                if self.option("group_repetition") == 'true':  # 选择组内重复时，样品数不能小于3
                    if df_count.iloc[i, 2] < 3:
                        raise OptionError('分组文件不满足参数要求')
                else:  # 选择组内不重复时，样品数为1
                    if df_count.iloc[i, 2] != 1:
                        raise OptionError('分组文件不满足参数要求')

    def line_chart_progress(self):
        if self.option('ranks') == 'row':  # 以行为标签时，表格不进行转置
            data_file = self.option('data_file').prop['path']
        else:
            data_file = self.output_dir + "/data_file_t.txt"
        f_data = pd.read_csv(data_file, sep='\t')
        f_group = pd.read_csv(self.option("group_file").prop["path"], sep='\t')
        f_data.rename(columns={f_data.columns[0]: 'sample_name'}, inplace=True)
        f_group.rename(columns={f_group.columns[0]: 'sample_name'}, inplace=True)
        if self.option('line_type') == "single":  # 选择为单折线时
            if self.option('group_repetition') == 'true':  # 选择进行组内重复
                df = pd.merge(f_group, f_data).drop('sample_name',axis=1)
                df_mean = df.iloc[:, 1].groupby(df.iloc[:, 0]).mean().reset_index(drop=False)
                df_mean.rename(columns={df_mean.columns[1]: 'mean'}, inplace=True)
                df_std =  df.iloc[:, 1].groupby(df.iloc[:, 0]).std().reset_index(drop=False)
                df_std.rename(columns={df_std.columns[1]: 'SD'}, inplace=True)
                df_out = pd.merge(df_mean, df_std)
                df_out.to_csv(os.path.join(self.output_dir, "group_result.txt"), header=True, index=0,
                              sep='\t', encoding="utf-8")
            else:  # 组内不重复
                df_out = pd.merge(f_group, f_data).drop('sample_name',axis=1)
                df_out.to_csv(os.path.join(self.output_dir, "group_result.txt"), header=True, index=0,
                              sep='\t', encoding="utf-8")
        else:  # 选择多折线
            if self.option('group_repetition') == 'true':  # 选择进行组内重复
                df = pd.merge(f_group, f_data).drop('sample_name',axis=1)
                df_mean = df.iloc[:, 2].groupby([df.iloc[:, 0], df.iloc[:, 1]]).mean().reset_index(drop=False)
                df_mean.rename(columns={df_mean.columns[2]: 'mean'}, inplace=True)
                df_std = df.iloc[:, 2].groupby([df.iloc[:, 0], df.iloc[:, 1]]).std().reset_index(drop=False)
                df_std.rename(columns={df_std.columns[2]: 'SD'}, inplace=True)
                df_out = pd.merge(df_mean, df_std)
                df_out.to_csv(os.path.join(self.output_dir, "group_result.txt"), header=True, index=0,
                              sep='\t', encoding="utf-8")
            else:
                df_out = pd.merge(f_group, f_data).drop('sample_name',axis=1)
                df_out.to_csv(os.path.join(self.output_dir, "group_result.txt"), header=True, index=0,
                              sep='\t', encoding="utf-8")

    def run(self):
        super(LineChartTool, self).run()
        self.file_check()
        self.line_chart_progress()
        self.end()
