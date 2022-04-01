# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20180911
# last modified: guhaidong

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd
from pymongo import MongoClient
from biocluster.config import Config
import shutil
import commands


class SelectTableAgent(Agent):
    """
    根据样本名和样本分组计算组和，
    """

    def __init__(self, parent):
        super(SelectTableAgent, self).__init__(parent)
        options = [
            {"name": "origin_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "select_table", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "group_method", "type": "string", "default": "average"},
            {"name": "scale", "type": "bool", "default": False},  ## 是否标准化,代谢用
            {"name": "select_origin_abu", "type": "outfile", "format": "sequence.profile_table"}, ## 有scale时使用
            {"name": "scale_method","type":"string","default":"UV"} #zouuguanqing 20190605
        ]
        self.add_option(options)
        self._memory_increase_step = 30  # modified by GHD @ 20180628


    def check_options(self):
        if not self.option("origin_table").is_set:
            raise OptionError("必须设置原始文件", code="34002501")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '16G'  #  by guhaidong @ 20180628
        # memory = 8 + 10 * self._rerun_time  # 每次重运行增加5G内存 by guhaidong @ 20180417
        # self._memory = "%sG" % memory

    def end(self):
        super(SelectTableAgent, self).end()


class SelectTableTool(Tool):
    def __init__(self, config):
        super(SelectTableTool, self).__init__(config)
        self.logger.info("按分组开始处理数据")
        self.abund_table_path = ""
        self.abund_table_percent_path = ""

    def cat_samples_percent(self):
        table = pd.DataFrame(pd.read_table(self.option("origin_table").prop["path"], sep='\t', index_col=0))
        table.columns = [str(i) for i in table.columns]  # guanqing.zou 20180512
        table_name = table.index.name
        table = table.ix[list((table > 0).any(axis=1))]  # 去除都为0的物种/功能/基因
        # if self.option("group_table").is_set:
        # if os.path.isfile(self.option("group_table").prop["path"]):  # by houshuang 20191010 is_set报错('bool' object has no attribute 'is_set')
        self.logger.info("group_table:{}".format(self.option("group")))
        if self.option("group").is_set:
            group = pd.DataFrame(pd.read_table(self.option("group").prop["path"], sep='\t', dtype={
                "#sample": "object"}))  # fix by qingchen.zhang 20191113 防止样本名称以数字的0开头
            # group["sample"] = group.index
            group = group.set_index("#sample")  # addqingchen.zhang 20191113 防止样本名称以数字的0开头
            group["sample"] = [str(i) for i in group.index]  # guanqing.zou 20180512
            group_sample = ""
            if self.option("group_method") == "sum":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).sum()  # 求和
            elif self.option("group_method") == "average":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).mean()  # 求均值
            elif self.option("group_method") == "middle":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).median()  # 中位数
            elif self.option("group_method") not in ["average", "sum", "middle"]:
                group_sample = group.join(table.T, on="sample")
                group_sample.drop(group_sample.columns[:2], axis=1, inplace=True)
            abund = group_sample.T
            abund.index.name = table_name
        else:
            abund = table
        abund['Col_sum'] = abund.apply(lambda x: x.sum(), axis=1)
        abund_table = abund.sort_values(by=['Col_sum'], ascending=0)
        del abund_table["Col_sum"]
        if len(abund_table) < 1:
            self.set_error('在所选参数下数据为空，请重新设置水平或分组方案参数!', code="32705401")
            self.set_error('在所选参数下数据为空，请重新设置水平或分组方案参数!', code="32705404")
        self.abund_table_path = os.path.join(self.output_dir, "metabolome.table.xls")
        new_abund_table = abund_table
        new_abund_table.to_csv(self.abund_table_path, sep="\t", encoding="utf-8")

        abund_table.columns = [str(i) for i in abund_table.columns]  ###guanqing.zou 20180514
        if  self.option("scale"):
            abund_table.loc['Row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
            sample_empt = []
            b = abund_table.apply(lambda x: x.sum(), axis=0)
            for i in range(len(abund_table.columns)):
                if b[i] > 0:
                    pass
                else:
                    sample_name = abund_table.columns[i]
                    sample_empt.append(sample_name)
            abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['Row_sum'], axis=1).drop('Row_sum')
            self.abund_table_percent_path = os.path.join(self.output_dir, "metabolome.percents.table.xls")
            new_abund_table_percent = abund_table_percent
            new_abund_table_percent.to_csv(self.abund_table_percent_path, sep="\t", encoding="utf-8")

    def run(self):
        super(SelectTableTool, self).run()
        self.cat_samples_percent()
        self.option("select_table").set_path(self.abund_table_path)
        if self.option("scale"):
            self.option("select_origin_abu").set_path(self.abund_table_percent_path)
        self.end()

