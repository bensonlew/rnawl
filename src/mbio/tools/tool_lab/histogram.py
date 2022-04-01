# !usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import os
import subprocess
import re


class HistogramAgent(Agent):
    """
    频率直方图
    """
    def __init__(self, parent):
        super(HistogramAgent, self).__init__(parent)
        options = [
            {"name": "input_table", "type": "infile", "format": "tool_lab.histogram_table"},
            {"name": "fq", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("input_table"):
            raise OptionError('必须设置作图数据')
        if not self.option("fq"):
            raise OptionError('必须设置频率')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "频率直方图数据"],
            ["./histogram.csv", "csv", "频率直方图数据"]
        ])
        super(HistogramAgent, self).end()


class HistogramTool(Tool):
    """
    频率直方图
    """
    def __init__(self, config):
        super(HistogramTool, self).__init__(config)
        self._version = 'v2.1-20140214'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def run(self):
        """
        运行
        """
        super(HistogramTool, self).run()
        self.get_table()
        self.set_output()
        # self.set_db()
        self.end()

    def get_table(self):
        """
        计算频率
        """
        n = 0
        fq_list = []
        fq = int(self.option("fq"))
        table = self.option("input_table").prop['path']
        with open(table, 'r') as f:
            list0 = f.readline().strip().split('\t')
            new_list = [eval(x) for x in list0]
            new_list.sort()
        while (n - 1) * fq < new_list[(len(new_list) - 1)]:
            fq_list.append(n * fq)
            n += 1
        list_cut = pd.cut(new_list, fq_list)
        new_list_cut = []
        for i in list_cut:
            if i not in new_list_cut:
                new_list_cut.append(i)
            else:
                pass
        df = pd.value_counts(list_cut)
        new_df = df[new_list_cut]
        new_df.to_csv(self.work_dir + "/new_table.csv", header=0, sep='\t')

    def set_output(self):
        """
        设置输出文件路径
        """
        path1 = self.output_dir + "/histogram.csv"
        if os.path.exists(path1):
            os.remove(path1)
        os.link(self.work_dir + "/new_table.csv", self.output_dir + "/histogram.csv")

    def set_db(self):
        """
        保存结果到Mongo库
        """
        self.logger.info("保存结果到Mongo")
        main_id = self.option("main_id")
        histogram_api = self.api.api("tool_lab.histogram_api")
        histogram_api.add_histogram_detail(main_id, self.output_dir)
        self.logger.info("保存结果到Mongo结束")
