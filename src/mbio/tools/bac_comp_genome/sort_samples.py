# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import linecache
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SortSamplesAgent(Agent):
    """
    传入一个group表，以及是否进行样本合并的参数生成一张丰度表并对并依照group表OTU表进行筛选合并
    """
    def __init__(self, parent):
        super(SortSamplesAgent, self).__init__(parent)
        options = [
            {"name": "abundance_file", "type": "infile", "format": "meta.otu.otu_table"},  # 输入层级的文件带有功能水平的
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # 输入的group表
            {"name": "method", "type": "string", "default": ""},  # 样本的合并方式, ""为不进行合并
            {"name": "out_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的结果OTU表
            {"name": "persent_otu_table", "type": "outfile", "format": "meta.otu.otu_table"},  # 输出的结果OTU表(百分比）
            {"name": "sample_del_warn", "type": "string", "default": "F"},  # 当样品的数据全为0时，是否报错提示
            {"name": "variable_del", "type": "string", "default": "T"},  #  是否去除都为0的物种/功能/基因
            {"name": "top", "type": "int", "default": ""},  # 是否取top物种/功能
        ]
        self.add_option(options)
        self.step.add_steps("sort_samples")
        self.on('start', self.start_sort_samples)
        self.on('end', self.end_sort_samples)

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
        if not self.option("abundance_file").is_set:
            raise OptionError("输入的丰度文件不能为空")
        # if not self.option("anno_file").is_set:
        #     raise OptionError("输入的注释文件不能为空")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["taxa.table.xls", "xls", "各样本物种丰度结果表"]
        ])
        super(SortSamplesAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "20G"


class SortSamplesTool(Tool):
    def __init__(self, config):
        super(SortSamplesTool, self).__init__(config)
        self.logger.info("SortSamplesMg按分组开始处理数据")
        self.abund_table_path = ""
        self.abund_table_percent_path = ""
        self.convert_dict = {
            "COG": "COG ID",
            "COG Function": "",
            "KO(KEGG Othology)": "",
            "KEGG Pathway Level3": "",
        }

    def get_level_abundance(self):
        """
        根据注释表和丰度表合并各层级的表
        :return:
        """
        abund_table = pd.read_table(self.option("abundance_file").prop["path"], sep='\t', index_col=0)
        # try:
        #     abund_table = abund_table.reset_index()
        # except:
        #     abund_table = abund_table
        abund_table.columns = [str(i) for i in abund_table.columns]  #防止样本名称为数字
        # abund_table.columns = ["sample"] + list(abund_table.columns)[1:]
        abund_table.index.name = "sample"
        table_name = abund_table.index.name
        table = abund_table.ix[list((abund_table > 0).any(axis=1))]  # 去除都为0的功能
        if self.option("group_table").is_set:
            group = pd.read_table(self.option("group_table").prop["path"], sep='\t', header=0)
            try:
                group = group.set_index("#sample")
            except:
                group = group.set_index("sample")
            group["sample"] = [str(i) for i in group.index]
            group_sample = ""
            if self.option("method") == "sum":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).sum()  # 求和
            elif self.option("method") == "average":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).mean()  # 求均值
                self.logger.info("group_sample: %s"%(group_sample.head()))
            elif self.option("method") == "middle":
                group_sample = group.join(table.T, on="sample").groupby(group.columns[0]).median()  # 中位数
            elif self.option("method") not in ["average", "sum", "middle"]:
                try:
                    all_group_list = list(group["#sample"])
                except:
                    all_group_list = list(group["sample"])
                # group_sample = table[[table.columns[0]] + [ x for x in table.columns[1:] if x in all_group_list]]
                group_sample = table[[ x for x in table.columns if x in all_group_list]]
            if self.option("method") != "":
                abund = group_sample.T
                abund.index.name = table_name
            else:
                abund = group_sample
        else:
            abund = table
        self.logger.info("abund: %s"%(abund.head()))
        abund['Col_sum'] = abund.apply(lambda x: x.sum(), axis=1)
        abund_table = abund.sort_values(by=['Col_sum'], ascending=0)
        del abund_table["Col_sum"]
        if len(abund_table) < 1:
            self.set_error('在所选参数下数据为空，请重新设置水平或分组方案参数!')

        self.abund_table_path = os.path.join(self.work_dir, "taxa.table.xls")
        if self.option("top") != "":
            top = int(self.option("top"))
            new_abund_table = abund_table.head(top)
        else:
            new_abund_table = abund_table
        new_abund_table.to_csv(self.abund_table_path, sep="\t", encoding="utf-8")#获取按照分组情况合并的绝对丰度表

        abund_table.loc['Row_sum'] = abund_table.apply(lambda x: x.sum(), axis=0)
        abund_table_percent = abund_table.apply(lambda x: x / abund_table.loc['Row_sum'], axis=1).drop('Row_sum')
        self.abund_table_percent_path = os.path.join(self.work_dir, "taxa.percents.table.xls")
        if self.option("top") != "":
            top = int(self.option("top"))
            new_abund_table_percent = abund_table_percent.head(top)
        else:
            new_abund_table_percent = abund_table_percent
        new_abund_table_percent.to_csv(self.abund_table_percent_path, sep="\t", encoding="utf-8")# 获得按照分组情况合并的百分比丰度表

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("正在生成结果文件目录")
        if os.path.exists(self.output_dir + "/taxa.table.xls"):
            os.remove(self.output_dir + "/taxa.table.xls")
        os.link(self.work_dir + "/taxa.table.xls", self.output_dir + "/taxa.table.xls")
        self.option("out_otu_table", self.output_dir + "/taxa.table.xls")
        if os.path.exists(self.output_dir + "/taxa.percents.table.xls"):
            os.remove(self.output_dir + "/taxa.percents.table.xls")
        os.link(self.work_dir + "/taxa.percents.table.xls", self.output_dir + "/taxa.percents.table.xls")
        self.option("persent_otu_table", self.output_dir + "/taxa.table.xls")


    def run(self):
        super(SortSamplesTool, self).run()
        self.get_level_abundance()
        self.set_output()
        self.end()
