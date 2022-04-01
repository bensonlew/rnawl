# -*- coding: utf-8 -*-
# __author__ = zengjing
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import re
import os


class CircosAgent(Agent):
    """
    小工具弦图：对二维表格(行为样本)进行合并
    """
    def __init__(self, parent):
        super(CircosAgent, self).__init__(parent)
        options = [
            {"name": "data_table", "type": "infile", "format": "toolapps.table"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},
            {"name": "combined_value", "type": "string", "default": "0%"}  #合并小于此数值的区域的值
        ]
        self.add_option(options)
        self.step.add_steps("circos")

    def check_options(self):
        if not self.option("data_table").is_set:
            raise OptionError("缺少输入的数据表格")
        if self.option("group_table").is_set:
            sample_names = self.option('group_table').prop['sample_name']
            for i in sample_names:
                if i not in self.option('data_table').prop['col_sample']:
                    raise OptionError('分组文件中的样本{}不存在于表格中，查看是否是数据表行列颠倒'.format(i))

    def set_resource(self):
        self._cpu = 3
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "Circos结果目录"],
        #     ["./result_data", "xls", "结果表"],
        # ])
        super(CircosAgent, self).end()


class CircosTool(Tool):
    def __init__(self, config):
        super(CircosTool, self).__init__(config)
        self._version = 1.0

    def get_group_detail(self):
        """根据分组文件得到具体的分组方案"""
        group_samples = {}  # 分组对应中新样本对应的旧样本
        with open(self.option('group_table').prop['new_table'], "r") as f:
            line = f.readline().rstrip()
            line = re.split("\t", line)
            for i in range(1, len(line)):
                group_samples[line[i]] = {}
            for item in f:
                item = item.rstrip().split("\t")
                for i in range(1, len(line)):
                    try:
                        if item[i] and item[i] not in group_samples[line[i]]:
                            group_samples[line[i]][item[i]] = []
                            group_samples[line[i]][item[i]].append(item[0])
                        elif item[i]:
                            group_samples[line[i]][item[i]].append(item[0])
                        else:
                            self.set_error("{}样本不在分组方案{}内".format(item[0], line[i]))
                    except:
                        self.set_error("{}样本不在分组方案{}内".format(item[0], line[i]))
        if not group_samples:
            self.set_error("分组方案文件格式错误，请检查")
        return group_samples

    def get_group_data_table(self):
        with open(self.option("data_table").prop["new_table"], "r") as f:
            lines = f.readlines()
            line = lines[0].rstrip().split("\t")
            if self.option('group_table').is_set:
                group_samples = self.get_group_detail()
                self.logger.info(group_samples)
                for group in group_samples:
                    new_samples = list(group_samples[group].keys())
                    new_sample_index = {}
                    for s in new_samples:
                        old_samples = group_samples[group][s]
                        new_sample_index[s] = []
                        for o in old_samples:
                            for i in range(len(line)):
                                if o == line[i]:
                                    new_sample_index[s].append(i)
                    table_path = os.path.join(self.work_dir, group + "_table.xls")
                    with open(table_path, "w") as w:
                        header = line[0] + "\t" + '\t'.join(new_samples) + "\n"
                        w.write(header)
                        for item in lines[1:]:
                            item = item.rstrip().split("\t")
                            w.write(item[0] + "\t")
                            tmp = []
                            for s in new_samples:
                                summary = 0
                                for i in new_sample_index[s]:
                                    summary += float(item[i])
                                tmp.append(str(summary))
                            w.write('\t'.join(tmp) + "\n")
                    self.combined_value_stat(fp=table_path, combined_value=self.option("combined_value"), out=group + "_result_data")
            else:
                table_path = os.path.join(self.work_dir, "all_table.xls")
                with open(table_path, "w") as w:
                    w.write(line[0] + "\t" + "all" + "\n")
                    for item in lines[1:]:
                        item = item.rstrip().split("\t")
                        w.write(item[0] + "\t")
                        summary = 0
                        for s in item[1:]:
                            summary += float(s)
                        w.write(str(summary) + "\n")
                self.combined_value_stat(fp=table_path, combined_value=self.option("combined_value"), out="all_result_data")

    def combined_value_stat(self, fp, combined_value, out):
        """对小于此数值的区域进行合并"""
        if combined_value.endswith("%"):
            combined_value = float(combined_value.split("%")[0]) / 100
        data = pd.read_table(fp, header=0)
        col_sum = data.sum(axis=1)
        col_value = col_sum.values
        row_sum = data.T.iloc[1:].sum(axis=1)
        row_value = row_sum.values
        name = data.T.iloc[0].to_frame(name="name")
        with open(fp, "r") as f, open(out, "w") as w:
            lines = f.readlines()
            w.write(lines[0])
            others = {}
            samples = lines[0].strip().split("\t")
            for s in samples[1:]:
                others[s] = 0
            other = False
            for line in lines[1:]:
                line = line.strip().split("\t")
                flag = False
                for i in range(1, len(line)):
                    percent = float(line[i]) / float(row_value[i-1])
                    if percent > float(combined_value):
                        flag = True
                if flag:
                    w.write('\t'.join(line) + "\n")
                else:
                    other = True
                    for i in range(1, len(line)):
                        for j in range(1, len(samples)):
                            if i == j:
                                others[samples[j]] += float(line[i])
            if other:
                other_value = []
                for s in samples[1:]:
                    other_value.append(str(others[s]))
                w.write("others" + "\t" + '\t'.join(other_value) + "\n")

    def set_output(self):
        file_names = os.listdir(self.work_dir)
        for name in file_names:
            if name.endswith("result_data"):
                os.link(os.path.join(self.work_dir, name), os.path.join(self.output_dir, name))

    def run(self):
        super(CircosTool, self).run()
        self.get_group_data_table()
        self.set_output()
        self.end()
