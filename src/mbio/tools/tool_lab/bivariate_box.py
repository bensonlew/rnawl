# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'yifei.fang'@ 20210824
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class BivariateBoxAgent(Agent):
    """
    BivariateBoxAgent
    """
    def __init__(self, parent):
        super(BivariateBoxAgent, self).__init__(parent)
        options = [
            {"name": "data_table", "type": "infile", "format": "tool_lab.table"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "sample_label", "type": "string", "default": "row"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("data_table").is_set:
            raise OptionError("必须输入table文件")
        if not self.option("group_table").is_set:
            raise OptionError("必须输入group文件")


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
        super(BivariateBoxAgent, self).end()


class BivariateBoxTool(Tool):
    """
    箱线数据计算
    """
    def __init__(self, config):
        super(BivariateBoxTool, self).__init__(config)

    def formattable(self, tablepath):
        """
        转置表格
        """
        this_table = tablepath
        if self.option('sample_label') != 'row':
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

    def calculate_bivariate_box(self):
        """
        开始处理数据表
        """
        f1 = pd.read_table(self.option("group_table").prop["path"], sep="\t")
        f2 = pd.read_table(self.formattable(self.option("data_table").prop["path"]), sep="\t")
        f1.columns = ["name", "group1", "group2"]
        f2.columns = ["name", "value"]
        data = pd.merge(f1, f2, on="name")
        group_data = data.groupby(["group1", "group2"])

        with open(self.option("group_table").prop["path"], "r") as f:
            lines = f.readlines()
            line = lines[0].strip().split("\t")

        with open(os.path.join(self.output_dir, "box_file.txt"), "w") as w:
            w.write(line[1] + "\t" + line[2] + "\t" + "min" + "\t" + "q1" + "\t" + "median" + "\t" + "q3" + "\t" + "max" + "\t" + "outlier" + "\n")
            error_dot_name_list = []
            for key, value in group_data:
                group1 = list(set(value["group1"].tolist()))[0]
                group2 = list(set(value["group2"].tolist()))[0]

                res = value["value"].describe()
                iqr = res['75%'] - res['25%']
                up_out_line = res['75%'] + iqr * 1.5
                down_out_line = res['25%'] - iqr * 1.5

                error_dot = value["value"][value["value"].map(lambda x: True if x > up_out_line or x < down_out_line else False)]
                error_dot1 = list(error_dot)
                list_error = []
                for i in error_dot1:
                    j = str(i)
                    list_error.append(j)

                if len(list_error) != 0:
                    error_dot_name = error_dot.index[0]
                    if error_dot_name not in error_dot_name_list:
                        error_dot_name_list.append(error_dot_name)
                else:
                    pass

                normal_dot = value["value"][value["value"].map(lambda x: False if x > up_out_line or x < down_out_line else True)]
                min = normal_dot.min()
                max = normal_dot.max()
                if min >= res['25%']:
                    nor_min = ""
                else:
                    nor_min = min
                if max <= res['75%']:
                    nor_max = ""
                else:
                    nor_max = max
                w.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(group1, group2, nor_min, res['25%'], res['50%'], res['75%'], nor_max, ",".join(list_error)))

            scatter_error = data.T.iloc[:, error_dot_name_list]
            scatter_error.T.to_csv(os.path.join(self.output_dir, "scatter_file.txt"), sep='\t', index=False)

    def run(self):
        """
        运行
        """
        super(BivariateBoxTool, self).run()
        self.calculate_bivariate_box()
        self.end()

