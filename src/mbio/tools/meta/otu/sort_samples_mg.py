# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

import os
import linecache
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.statistical.normalization import Normalization


class SortSamplesMgAgent(Agent):
    """
    传入一个group表，以及是否进行样本合并的参数生成一张丰度表并对并依照group表OTU表进行筛选合并
    """
    def __init__(self, parent):
        super(SortSamplesMgAgent, self).__init__(parent)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入丰度文件
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # 输入的group表
            {"name": "method", "type": "string", "default": ""},  # 样本的合并方式, ""为不进行合并
            {"name": "out_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的结果OTU表
            {"name": "level_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的结果OTU表(百分比）
            {"name": "others", "type": "float", "default": ""},  # 组成分析中用于将丰度小于0.01/其它的物种归为others
            {"name": "top", "type": "int", "default": ""},  # 热图取top物种/功能
            {"name": "sample_del_warn", "type": "string", "default": "F"},  # 当样品的数据全为0时，是否报错提示
            {"name": "variable_del", "type": "string", "default": "T"},  #  是否去除都为0的物种/功能/基因
            {"name": "get_only_percent", "type" : "string", "default" : "F"},  # T 时 计算相对丰度，不做other处理。注意默认值不能改成T
            {"name": "normalization", "type": "string", "default": ""},  # 对数据进行标准化, row, col, 或 “”
        ]
        self.add_option(options)
        self.step.add_steps("sort_samples")
        self.on('start', self.start_sort_samples)
        self.on('end', self.end_sort_samples)
        self._memory_increase_step = 40  # 每次重运行增加内存40G by qingchen.zhang @ 20190916

    def start_sort_samples(self):
        self.step.sort_samples.start()
        self.step.update()

    def end_sort_samples(self):
        self.step.sort_samples.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option("in_otu_table").is_set:
            raise OptionError("输入的丰度文件不能为空", code="32705401")
        if self.option("method"):
            if self.option("method") not in ["", "no", "none", "No", "None", None, "average", "sum", "middle"]:
                raise OptionError("参数method设置错误！", code="32705402")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("others") == "":
            result_dir.add_relpath_rules([
                [".", "", "结果输出目录"],
                ["taxa.table.xls", "xls", "各样本物种丰度结果表"]
            ])
        else:
            result_dir.add_relpath_rules([
                [".", "", "结果输出目录"],
                ["taxa.table.xls", "xls", "各样本物种丰度结果表"],
                ["taxa.percents.table.xls", "xls", "各样本物种相对丰度结果表"]
            ])
        super(SortSamplesMgAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "10G"


class SortSamplesMgTool(Tool):
    def __init__(self, config):
        super(SortSamplesMgTool, self).__init__(config)
        self.logger.info("SortSamplesMg按分组开始处理数据")
        self.abund_table_path = ""
        self.abund_table_percent_path = ""

    def cat_samples_percent(self):
        table = pd.DataFrame(pd.read_table(self.option("in_otu_table").prop["path"], sep='\t', index_col=0))
        table.columns = [str(i) for i in table.columns] #guanqing.zou 20180512
        table_name = table.index.name
        table = table.ix[list((table > 0).any(axis=1))]  # 去除都为0的物种/功能/基因
        # if self.option("group_table").is_set:
        # if os.path.isfile(self.option("group_table").prop["path"]):  # by houshuang 20191010 is_set报错('bool' object has no attribute 'is_set')
        self.logger.info("group_table:{}".format(self.option("group_table")))
        if self.option("group_table"):
            group = pd.DataFrame(pd.read_table(self.option("group_table").prop["path"], sep='\t',dtype={"#sample":"object"})) # fix by qingchen.zhang 20191113 防止样本名称以数字的0开头
            #group["sample"] = group.index
            group = group.set_index("#sample") #addqingchen.zhang 20191113 防止样本名称以数字的0开头
            group["sample"] = [str(i) for i in group.index] #guanqing.zou 20180512
            group_sample = ""
            if self.option("method") == "sum":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).sum()  # 求和
            elif self.option("method") == "average":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).mean()  # 求均值
            elif self.option("method") == "middle":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).median()  # 中位数
            elif self.option("method") not in ["average", "sum", "middle"]:
                group_sample = group.join(table.T, on="sample")
                group_sample.drop(group_sample.columns[:2], axis=1, inplace=True)
            abund = group_sample.T
            abund.index.name = table_name
        else:
            abund = table
        abund['Col_sum'] = abund.apply(lambda x: x.sum(), axis=1)
        abund_table = abund.sort_values(by=['Col_sum'], ascending=0)
        del abund_table["Col_sum"]
        if self.option("variable_del") == "T":
            abund_table = abund_table.ix[list((abund_table > 0).any(axis=1))]  # 去除都为0的物种/功能/基因
        elif self.option("variable_del") == "F":
            abund_table = abund_table
        if len(abund_table) < 1:
            self.set_error('在所选参数下数据为空，请重新设置水平或分组方案参数!', code="32705401")
            self.set_error('在所选参数下数据为空，请重新设置水平或分组方案参数!', code="32705404")
        self.abund_table_path = os.path.join(self.output_dir, "taxa.table.xls")
        if self.option("top") != "":
            top = int(self.option("top"))
            new_abund_table = abund_table.head(top)
        else:
            new_abund_table = abund_table
        if self.option('normalization'):
            new_abund_table = Normalization(new_abund_table,
                                            by=self.option('normalization'),
                                            norm_meth='z').run()
        new_abund_table.to_csv(self.abund_table_path, sep="\t", encoding="utf-8")

        abund_table.columns = [str(i) for i in abund_table.columns] ###guanqing.zou 20180514
        if self.option("others") != "" or self.option("top") != "" or self.option("sample_del_warn") != "F" or self.option("get_only_percent") == 'T':
            abund_table.loc['Row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
            sample_empt = []
            b = abund_table.apply(lambda x: x.sum(), axis=0)
            for i in range(len(abund_table.columns)):
                if b[i] > 0:
                    pass
                else:
                    sample_name = abund_table.columns[i]
                    sample_empt.append(sample_name)
            if sample_empt:
                try:  # 增加判断，如果是工作流调用的此程序，不中断报错，自动跳过 by GHD @ 20180320
                    if self.get_workflow()._sheet.UPDATE_STATUS_API in ["tsanger", "sanger"]:
                        self.logger.info("工作流跳过此分析")
                        self.end()
                    else:
                        self.set_error('样品：%s在所选参数下数据均为0，请剔除该样品或重新设置参数!', variables=(sample_empt), code="32705402")
                except:
                    self.set_error('样品：%s在所选参数下数据均为0，请剔除该样品或重新设置参数!', variables=(sample_empt), code="32705403")
                # self.set_error('样品：%s在所选参数下数据均为0，请剔除该样品或重新设置参数!' % sample_empt)
            abund_table_percent = abund_table.apply(lambda x: x/abund_table.loc['Row_sum'], axis=1).drop('Row_sum')
            self.abund_table_percent_path = os.path.join(self.output_dir, "taxa.percents.table.xls")
            if self.option("top") != "": ## 已经排过顺序了
                top = int(self.option("top"))
                new_abund_table_percent = abund_table_percent.head(top)
            else:
                new_abund_table_percent = abund_table_percent
            new_abund_table_percent.to_csv(self.abund_table_percent_path, sep="\t", encoding="utf-8")

    def get_others(self, abund_table_percent_path):  # add by zhujuan 2017.10.12
        df = pd.DataFrame(pd.read_table(abund_table_percent_path, sep='\t', index_col=0))
        new_df = df.ix[list((df > self.option("others")).any(axis=1))]
        new_df2 = new_df.copy()
        others = df.ix[list((df < self.option("others")).all(axis=1))]
        if len(others) > 0:
            new_df2.loc["others"] = others.apply(lambda x: x.sum(), axis=0)
        other = os.path.join(self.output_dir, "taxa.percents.table.xls")
        new_df2.to_csv(other, sep="\t", encoding="utf-8")

    def run(self):
        super(SortSamplesMgTool, self).run()
        self.cat_samples_percent()
        self.option("out_otu_table").set_path(self.abund_table_path)
        #self.option("level_otu_table").set_path(self.abund_table_percent_path)
        if self.option("get_only_percent") == 'T':
            self.end()

        if self.option("others") != "":
            self.option("level_otu_table").set_path(self.abund_table_percent_path)
            self.get_others(self.abund_table_percent_path)
        self.end()
