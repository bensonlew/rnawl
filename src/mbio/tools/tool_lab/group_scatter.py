# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
import os
import numpy as np
import pandas as pd
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class GroupScatterAgent(Agent):
    def __init__(self, parent):
        super(GroupScatterAgent, self).__init__(parent)
        options = [
            {"name": "data_file", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group_file", "type": "infile", "format": "toolapps.group_table"},
            {"name": "ranks", "type": "string", "default": "row"},      #选择行、列标签
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        self.logger.info(self.option('data_file'))
        if not self.option("data_file").is_set:
            raise OptionError("必须设置输入数据文件")
        if not self.option("group_file").is_set:
            raise OptionError("必须设置输入分组文件")
        with open(self.option("group_file").path,'r') as f1:  #检查分组中的样本数  # 需要添加一个样本对比
            group_dict = {}
            for i in f1.readlines()[1:]:
                items = i.strip().split('\t')
                if items[1] not in group_dict:
                    group_dict[items[1]] = [items[0]]
                else:
                    group_dict[items[1]].append(items[0])
            for i in group_dict.keys():
                if len(group_dict[i]) < 4:
                    raise OptionError("每组中样本数不少于4")
        with open(self.option("data_file").path,"r") as f2:
            line = f2.readlines()
            if self.option('ranks') == 'row':
                if len(line[0].strip().split("\t")) != 2:
                    raise OptionError("行标签或列标签选择错误")
            else:
                if len(line) != 2:
                    raise OptionError("行标签或列标签选择错误")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(GroupScatterAgent, self).end()


class GroupScatterTool(Tool):
    def __init__(self, config):
        super(GroupScatterTool, self).__init__(config)

    def t_table(self, table_file, new_table):
        """
        转置表格实现
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)
        
    def group_scatter_progress(self):
        if self.option('ranks') == 'row':     #根据行列标签参数决定是否对表格进行转置
            data_table = self.option('data_file').prop['path']
        else:
            data_table = self.output_dir + "/data_file_t.txt"
            self.t_table(self.option('data_file').prop['path'] , data_table)
        with open(self.option("group_file").prop["path"]) as f1, \
             open(data_table) as f2, \
             open(os.path.join(self.output_dir, "group_result.txt"), "w+") as w:
            dict_group = {}
            dict_data = {}
            list_group = []
            for line in f1.readlines()[1:]:
                item = line.strip().split("\t")
                if item[1] not in dict_group:
                    list_group.append(item[1])
                    dict_group[item[1]] = [item[0]]
                else:
                    dict_group[item[1]].append(item[0])
            for lines in f2.readlines()[1:]:
                samples = lines.strip().split("\t")
                if samples[0] not in dict_data.keys():
                    dict_data[samples[0]] = samples[1]
                else:
                    pass
            self.logger.info(dict_data)
            self.logger.info(list_group)
            w.write("group_name\tmin\tq1\tmedian\tq3\tmax\r\n")  #计算输出每组内数据的最小值，1/4分位值，中位值，3/4分位值，已经最大值
            for group_name in list_group:
                list = [float(dict_data[x]) for x in dict_group[group_name]]
                w.write(group_name + "\t" + str(round(min(list),2)) + "\t" + str(round(np.percentile(list,25),2)) + "\t"\
                + str(round(np.median(list),2)) + "\t" + str(round(np.percentile(list,75),2)) + "\t" + str(round(max(list),2)) + "\r\n")

    def run(self):
        '''
        运行
        :return:
        '''
        super(GroupScatterTool, self).run()
        self.group_scatter_progress()
        self.end()
