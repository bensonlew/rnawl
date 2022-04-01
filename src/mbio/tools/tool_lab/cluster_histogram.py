# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd


class ClusterHistogramAgent(Agent):
    def __init__(self, parent):
        super(ClusterHistogramAgent, self).__init__(parent)
        options = [
            {"name": "data_file", "type": "infile", "format": "tool_lab.table"}, #输入数据表
            {"name": "group_file", "type": "infile", "format": "tool_lab.group_table"}, #输入分组文件
            {"name": "ranks", "type": "string", "default": "row"}, #行、列标签，默认为行
            {"name": "group_repetition","type":"string","default":"true"}, #组内重复，true or false
            {'name': 'error_bar', 'type': 'string', "default": "std"},  #标准差或标准误,std or sem
            ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        df = pd.read_csv(self.option("group_file").prop["path"], sep='\t')
        if len(df.columns) != 3:
            raise OptionError('分组文件必须为三列')
        df_count = df.iloc[:, 0].groupby([df.iloc[:, 1], df.iloc[:, 2]]).count().reset_index(drop=False)  #计算分组中样本数
        for i in df_count.index:
            if self.option("group_repetition") == 'true': # 组内重复
                if df_count.iloc[i,2] < 2:
                    raise OptionError('分组中样本数不符合要求，请检查')
            else: # 组内不重复
                if df_count.iloc[i,2] != 1:
                    raise OptionError('分组中样本数不符合要求，请检查')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(ClusterHistogramAgent, self).end()


class ClusterHistogramTool(Tool):
    def __init__(self, config):
        super(ClusterHistogramTool, self).__init__(config)

    def t_table(self, table_file, new_table):
        """
        转置表格
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)

    def data_process(self):
        """
        数据处理计算
        """
        if self.option('ranks') == 'row':  # 选择样本名以行为标签不进行转置
            data_table = self.option('data_file').prop['path']
        else:
            data_table = self.output_dir + "/data_file_t.txt"
            self.t_table(self.option('data_file').prop['path'] , data_table)
        f_data = pd.read_csv(data_table, sep='\t')
        f_group = pd.read_csv(self.option("group_file").prop["path"], sep='\t')
        f_data.rename(columns={f_data.columns[0]: 'sample_name'}, inplace=True)
        f_group.rename(columns={f_group.columns[0]: 'sample_name'}, inplace=True)
        if self.option('group_repetition') == 'true':  # 选择组内重复
            df = pd.merge(f_group, f_data).drop('sample_name', axis=1)
            df_mean = df.iloc[:, 2].groupby([df.iloc[:, 0], df.iloc[:, 1]]).mean().reset_index(drop=False)
            df_mean.rename(columns={df_mean.columns[2]: 'mean'}, inplace=True)
            df_mean.to_csv(os.path.join(self.output_dir, "merge_result.txt"), header=True, index=0, sep='\t', encoding="utf-8")
            if self.option("error_bar") == "std":  # 计算标准差std
                df_std = df.iloc[:, 2].groupby([df.iloc[:, 0], df.iloc[:, 1]]).std().reset_index(drop=False)
                df_std.rename(columns={df_std.columns[2]: 'SD'}, inplace=True)
                df_out = pd.merge(df_mean, df_std)
                df_out.to_csv(os.path.join(self.output_dir, "group_result.txt"),header=True,index=0,sep='\t', encoding="utf-8")
            else:  # 计算标准误sem
                df_sem = df.iloc[:, 2].groupby([df.iloc[:, 0], df.iloc[:, 1]]).sem().reset_index(drop=False)
                df_sem.rename(columns={df_sem.columns[2]: 'SE'}, inplace=True)
                df_out = pd.merge(df_mean, df_sem)
                df_out.to_csv(os.path.join(self.output_dir, "group_result.txt"),header=True,index=0,sep='\t', encoding="utf-8")
        else:  # 选择组内不重复
            df_out = pd.merge(f_group, f_data).drop('sample_name', axis=1)
            df_out.to_csv(os.path.join(self.output_dir, "merge_result.txt"), header=True, index=0, sep='\t', encoding="utf-8")

    def run(self):
        super(ClusterHistogramTool, self).run()
        self.data_process()
        self.end()