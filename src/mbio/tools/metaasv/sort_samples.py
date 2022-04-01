# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import linecache
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import filter_zero_and_replace


class SortSamplesAgent(Agent):
    """
    传入一个group表，以及是否进行样本合并的参数生成一张丰度表并对并依照group表OTU表进行筛选合并
    """
    def __init__(self, parent):
        super(SortSamplesAgent, self).__init__(parent)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入丰度文件
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # 输入的group表
            {"name": "method", "type": "string", "default": ""},  # 样本的合并方式, ""为不进行合并
            {"name": "out_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的结果OTU表
            {"name": "level_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的结果OTU表(百分比）
            {"name": "others", "type": "float", "default": ""},  # 组成分析中用于将丰度小于0.01/其它的物种归为others
            {"name": "top", "type": "int", "default": 0},  # 热图取top物种/功能
            {"name": "fill_zero", "type": "string", "default": "true"},  #  是否对0全表用最小值进行填充
            {"name": "abundance_method", "type" : "string", "default" : "all"}  #计算相对丰度和绝对丰度
        ]
        self.add_option(options)
        self.step.add_steps("sort_samples")
        self.on('start', self.start_sort_samples)
        self.on('end', self.end_sort_samples)
        self._memory_increase_step = 40

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
            raise OptionError("输入的丰度文件不能为空")
        if self.option("method"):
            if self.option("method") not in ["", "no", "none", "No", "None", None, "average", "sum", "middle"]:
                raise OptionError("参数method设置错误！")

    def end(self):
        super(SortSamplesAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "10G"


class SortSamplesTool(Tool):
    def __init__(self, config):
        super(SortSamplesTool, self).__init__(config)
        self.logger.info("SortSamplesMg按分组开始处理数据")
        self.abund_table_path = ""
        self.abund_table_percent_path = ""

    def cat_samples_percent(self):
        """
        功能：根据分组文件合并丰度表，计算百分比
        :return:
        """
        table = pd.DataFrame(pd.read_table(self.option("in_otu_table").prop["path"], sep='\t', index_col=0, header=0))
        table.columns = [str(i) for i in table.columns]
        table_name = table.index.name
        if self.option("group_table"):
            group_table = pd.DataFrame(pd.read_table(self.option("group_table").prop["path"], sep='\t',dtype={"#sample":"object"}))
            group_table.columns = [str(i) for i in group_table.columns]
            if len(group_table.columns) >2:
                new_columns = list(group_table.columns)[0:2]
                group_table = group_table[new_columns]
            try:
                group_table = group_table.set_index("#sample")
            except:
                columns_list = group_table.columns
                group_index_name =  columns_list[0]
                group_table = group_table.set_index(group_index_name)
            group_table["sample"] = [str(i) for i in group_table.index]
            group_sample = ""
            if self.option("method") in ["sum"]:
                group_sample = group_table.join(table.T, on="sample").groupby(group_table.columns[0]).sum()  # 求和
            elif self.option("method") in ["average"]:
                group_sample = group_table.join(table.T, on="sample").groupby(group_table.columns[0]).mean()  # 求均值
            elif self.option("method") in ["middle"]:
                group_sample = group_table.join(table.T, on="sample").groupby(group_table.columns[0]).median()  # 中位数
            elif self.option("method") not in ["average", "sum", "middle"]:
                group_sample = group_table.join(table.T, on="sample")
                group_sample.drop(group_sample.columns[:2], axis=1, inplace=True)
            abund = group_sample.T
            try:
                abund.drop(index="sample",axis=0, inplace=True)
            except:
                pass
            abund.index.name = table_name
            self.logger.info("工作流跳过此分析{}".format(abund.head()))
        else:
            abund = table
        abund['all_sum'] = abund.apply(lambda x: x.sum(), axis=1)
        abund_table = abund.sort_values(by=['all_sum'], ascending=0)
        del abund_table["all_sum"]
        abund_table = abund_table.ix[list((abund_table > 0).any(axis=1))]  # 去除都为0的物种
        if len(abund_table) < 1:
            raise Exception('在所选参数下数据为空，请重新设置水平或分组方案参数!')
        sample_empt = []
        self.logger.info("工作流跳过此分析{}".format(abund_table.head()))
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
                    self.set_error('样品：%s在所选参数下数据均为0，请剔除该样品或重新设置参数!')
            except:
                self.set_error('样品：%s在所选参数下数据均为0，请剔除该样品或重新设置参数!')

        dir_name = self.work_dir
        abundance_method = self.option("abundance_method")
        if abundance_method in ['absolute']:##只有绝对丰度
            abund_table_path = os.path.join(dir_name, "taxa.table.xls")
            if self.option("top") not in ["", 0]:##只有相对丰度
                top_num = int(self.option("top"))
                abund_table = abund_table.head(top_num)
            abund_table.to_csv(abund_table_path, sep="\t", encoding="utf-8")

        elif abundance_method in ['relative']:##只有相对丰度
            abund_table.columns = [str(i) for i in abund_table.columns]
            abund_table.loc['row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
            abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['row_sum'], axis=1).drop('row_sum')
            abund_table.drop("row_sum")
            self.abund_table_percent_path = os.path.join(dir_name, "taxa.percents.table.xls")
            abund_table_percent_path = os.path.join(dir_name, "taxa.percents.mongo.xls")
            if self.option("fill_zero") in ["true"]:
                fill_zero_path = os.path.join(dir_name, "fill_zero.xls")
                abund_table_percent.to_csv(fill_zero_path, sep="\t", encoding="utf-8")
                abund_table_percent2 = filter_zero_and_replace(fill_zero_path,abund_table_percent_path)##结果文件用于导表
                # self.logger.info("abund_table: {}".format(abund_table_percent.head()))
            if self.option("top") not in ["", 0]:
                top_num = int(self.option("top"))
                abund_table = abund_table.head(top_num)
            # abund_table_list = list(abund_table.iloc[:, 0])
            abund_table_list = list(abund_table.index)
            abund_table_percent = abund_table_percent[(abund_table_percent.index).isin(abund_table_list)]
            if self.option("others") != "":
                new_df = abund_table_percent.ix[list((abund_table_percent > self.option("others")).any(axis=1))]
                new_df2 = new_df.copy()
                others = abund_table_percent.ix[list((abund_table_percent < self.option("others")).all(axis=1))]
                if len(others) > 0:
                    new_df2.loc["others"] = others.apply(lambda x: x.sum(), axis=0)
                abund_table_percent = new_df2
            abund_table_percent.to_csv(self.abund_table_percent_path, sep="\t", encoding="utf-8")##结果文件用于提供给用户

        elif abundance_method in ['all']:##含有绝对丰度和相对丰度
            abund_table.columns = [str(i) for i in abund_table.columns]
            # self.logger.info("abund_table.columns: {}".format(abund_table.columns))
            abund_table.loc['row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
            abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['row_sum'], axis=1).drop('row_sum')
            abund_table.drop('row_sum',inplace=True)
            # self.logger.info("abund_table: {}".format(abund_table.head()))
            self.abund_table_percent_path = os.path.join(dir_name, "taxa.percents.table.xls")
            abund_table_percent_path = os.path.join(dir_name, "taxa.percents.mongo.xls")
            abund_table_path = os.path.join(dir_name, "taxa.table.xls")
            if self.option("fill_zero") in ["true"]:
                fill_zero_path = os.path.join(dir_name, "fill_zero.xls")
                abund_table_percent.to_csv(fill_zero_path, sep="\t", encoding="utf-8")
                abund_table_percent2 = filter_zero_and_replace(fill_zero_path,abund_table_percent_path)##结果文件用于导表
                # self.logger.info("abund_table: {}".format(abund_table_percent.head()))
            if self.option("top") not in ["", 0]:
                top_num = int(self.option("top"))
                abund_table = abund_table.head(top_num)
            # abund_table_list = list(abund_table.iloc[:,0])
            abund_table_list = list(abund_table.index)
            self.logger.info("abund_table_list: {}".format(abund_table_list))
            self.logger.info("abund_table_percent: {}".format(abund_table_percent.head()))
            abund_table.to_csv(abund_table_path, sep="\t", encoding="utf-8")
            # abund_table_percent = abund_table_percent[abund_table_percent.iloc[:,0].isin(abund_table_list)]
            abund_table_percent = abund_table_percent[(abund_table_percent.index).isin(abund_table_list)]
            if self.option("others") not in ["", 0]:
                new_df = abund_table_percent.ix[list((abund_table_percent > self.option("others")).any(axis=1))]
                new_df2 = new_df.copy()
                others = abund_table_percent.ix[list((abund_table_percent < self.option("others")).all(axis=1))]
                if len(others) > 0:
                    new_df2.loc["others"] = others.apply(lambda x: x.sum(), axis=0)
                abund_table_percent = new_df2
            abund_table_percent.to_csv(self.abund_table_percent_path, sep="\t", encoding="utf-8")##结果文件用于提供给用户

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        if os.path.exists(os.path.join(self.output_dir, "taxa.table.xls")):
            os.remove(os.path.join(self.output_dir, "taxa.table.xls"))
        os.link(os.path.join(self.work_dir, "taxa.table.xls"), os.path.join(self.output_dir, "taxa.table.xls"))
        if os.path.exists(os.path.join(self.output_dir, "taxa.percents.table.xls")):
            os.remove(os.path.join(self.output_dir, "taxa.percents.table.xls"))
        os.link(os.path.join(self.work_dir, "taxa.percents.table.xls"), os.path.join(self.output_dir, "taxa.percents.table.xls"))
        self.option("level_otu_table").set_path(os.path.join(self.output_dir, "taxa.percents.table.xls"))
        self.option("out_otu_table").set_path(os.path.join(self.output_dir, "taxa.table.xls"))
        if self.option("fill_zero") in ["true"]:
            if os.path.exists(os.path.join(self.output_dir, "taxa.percents.mongo.xls")):
                os.remove(os.path.join(self.output_dir, "taxa.percents.mongo.xls"))
            os.link(os.path.join(self.work_dir, "taxa.percents.mongo.xls"), os.path.join(self.output_dir, "taxa.percents.mongo.xls"))

    def run(self):
        super(SortSamplesTool, self).run()
        self.cat_samples_percent()
        self.set_output()
        self.end()
