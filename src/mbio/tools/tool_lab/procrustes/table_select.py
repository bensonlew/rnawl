# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import math
import os
import glob
import unittest
import re


class TableSelectAgent(Agent):
    """
    PCoA
    """
    def __init__(self, parent):
        super(TableSelectAgent, self).__init__(parent)
        options = [
            {"name": "origin_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "select_table", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "group_method", "type": "string", "default": "none"},
            {"name": "log", 'type': 'string', 'default': 'log2'},
            {"name": 'group_names', 'type': 'string'},
            {'name': 'method', 'type': 'string'},
        ]
        self.add_option(options)
        self.step.add_steps('table_select')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.table_select.start()
        self.step.update()

    def step_end(self):
        self.step.table_select.finish()
        self.step.update()

    def check_options(self):
        if not self.option("origin_table").is_set:
            raise OptionError("必须设置原始文件", code="34002501")
        if not self.option('group').is_set:
            raise OptionError("必须设置输入分组文件")
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 2
        self._memory = "16G"

    def end(self):
        super(TableSelectAgent, self).end()


class TableSelectTool(Tool):
    def __init__(self, config):
        super(TableSelectTool, self).__init__(config)
        self._version = "1.0"
        self.abund_table_path = os.path.join(self.output_dir, os.path.basename(self.option('origin_table').prop['path']))

    def dataframe_process(self):

        def log_transform(x, y):
            if x > 0:
                x = x
            else:
                x = 1
            out = math.log(x, y)
            return out

        group_dict = dict()
        samples = list()
        with open(self.option('group').prop['path'], 'r') as g:
            g.readline()
            for line in g:
                items = line.strip().split('\t')
                if items[1] not in group_dict.keys():
                    group_dict[items[1]] = list()
                group_dict[items[1]].append(items[0])
        if self.option('group_names'):
            groups = self.option('group_names').split(';')
        else:
            groups = group_dict.keys()
        for i in groups:
            samples += group_dict.get(i)
        df = pd.read_table(self.option("origin_table").prop["path"], sep='\t', index_col=0)
        try:
            df = df[samples]
        except Exception as e:
            diff = [i not in df.columns for i in samples]
            self.logger.info('Unable to find samples: {} in dataframe'.format(','.join(diff)))
            self.set_error('提取指定分组失败{}'.format(e))
        if self.option('log') == 'log2':
            df = df.applymap(lambda x: log_transform(x, 2))
        elif self.option('log') == 'log10':
            df = df.applymap(lambda x: log_transform(x, 10))
        return df

    def cat_samples_percent(self):
        table = pd.DataFrame(self.dataframe_process())
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
            if self.option('group_names'):
                groups = self.option('group_names').split(';')
                group = group[(group['group'].isin(groups))]
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
        new_abund_table = abund_table
        new_abund_table.to_csv(self.abund_table_path, sep="\t", encoding="utf-8")

    def quote_process(self):
        with open(self.abund_table_path, 'r') as r:
            raw = r.read()
        processed = raw.replace('\'', '')
        final = processed.replace('\"', '')
        with open(self.abund_table_path, 'w') as f:
            f.write(final)

    def run(self):
        super(TableSelectTool, self).run()
        if self.option('method') == 'pcoa':
            self.cat_samples_percent()
        else:
            new_abund_table = self.dataframe_process()
            new_abund_table.to_csv(self.abund_table_path, sep="\t", encoding="utf-8")
        self.quote_process()
        self.option("select_table").set_path(self.abund_table_path)
        self.end()

