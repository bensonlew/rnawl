# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last_modifiy:2017.10.09

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd


class CreateAbundTableLefseAgent(Agent):
    """
    生成各项分析的不同数据库的标准丰度表格
    """
    def __init__(self, parent):
        super(CreateAbundTableLefseAgent, self).__init__(parent)
        options = [
            {"name": "anno_table", "type": "infile", "format": "meta.profile"},  # 各数据库的注释表格
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            {"name": "level_type", "type": "string", "default": ""},
            # 注释表的字段，eg：Pathway，Level1，Level2
            {"name": "level_type_name", "type": "string", "default": ""},
            # 注释表字段的具体levelname，eg：Level1下Metabolism(对应的KO)
            {"name": "lowest_level", "type": "string", "default": ""},  # 注释表数据库对应的最低分类，eg：KEGG的ko
            {"name": "out_table", "type": "outfile", "format": "meta.otu.otu_table"},
        ]
        self.add_option(options)
        #self._memory_increase_step = 30  # 每次重运行增加30G内存 add by GHD @ 20181205
        self._memory_increase_step = 40  # 每次重运行增加内存40G by qingchen.zhang @ 20190916

    def check_options(self):
        if not self.option("anno_table").is_set and not self.option("gene_list").is_set:
            raise OptionError("请传入注释表格或者基因list文件！", code="32703901")
        if not self.option("geneset_table").is_set:
            raise OptionError("请传入基因丰度文件！", code="32703902")
        if self.option("level_type") == "" and self.option("anno_table").is_set:
            raise OptionError("请提供筛选的level水平", code="32703903")
        if self.option("level_type_name") != "":
            if self.option("lowest_level") == "":
                raise OptionError("请提供数据库对应的最低分类", code="32703904")

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G' # 内存5G增加到10G by GHD @20180205

    def end(self):
        super(CreateAbundTableLefseAgent, self).end()


class CreateAbundTableLefseTool(Tool):
    def __init__(self, config):
        super(CreateAbundTableLefseTool, self).__init__(config)

    def create_abund_table(self):
        geneset_table_path = self.option("geneset_table").prop["path"]
        self.logger.info("geneset_table is : " + geneset_table_path)
        geneset_table = pd.read_table(geneset_table_path, sep='\t', header=0)
        new_otu_file_path = os.path.join(self.output_dir, "new_abund_table.xls")
        if self.option("gene_list").is_set:
            gene_list_path = self.option("gene_list").prop["path"]
            self.logger.info("gene_list is : " + gene_list_path)
            gene_list = pd.read_table(gene_list_path, sep='\t', header=0)
            gene_table = geneset_table.set_index("GeneID").ix[list(gene_list["GeneID"])]
            gene_table.to_csv(new_otu_file_path, sep="\t")
        else:
            anno_table_path = self.option("anno_table").prop["path"]
            self.logger.info("anno table is : " + anno_table_path)
            anno_table = pd.read_table(anno_table_path, sep='\t', header=0)
            a = pd.DataFrame(anno_table)
            if self.option("level_type_name") != "":
                a = a.ix[:, ["#Query", self.option("lowest_level"), self.option("level_type")]]
                a.columns = ["GeneID", self.option("lowest_level"), self.option("level_type")]
            else:
                a = a.ix[:, ["#Query", self.option("level_type")]]
                a.columns = ["GeneID", self.option("level_type")]
            b = pd.DataFrame(geneset_table)
            abund = a.merge(b, on='GeneID', how='inner')
            abund[self.option('level_type')] = abund[self.option('level_type')].apply(lambda x: x.replace('; ', ';'))
            level_type_abund = abund.drop(self.option("level_type"), axis=1).join(
                abund[self.option("level_type")].str.split(';', expand=True).stack().reset_index(
                    level=1, drop=True).rename(self.option("level_type")))
            if self.option("level_type_name") != "":
                if self.option("lowest_level") in ["Pathway"]:
                    abund[self.option('lowest_level')] = abund[self.option('lowest_level')].apply(lambda x: x.replace('; ', ';'))
                    pathway = abund[self.option("lowest_level")].str.split(";", expand=True).stack().reset_index(
                        level=1, drop=True).rename(self.option("lowest_level"))
                    level_type_abund[self.option("lowest_level")] = list(pathway)
                    level_type_abund = level_type_abund[level_type_abund[self.option("level_type")] == self.option(
                        "level_type_name")]
                    level_type_abund = level_type_abund[level_type_abund[self.option("lowest_level")] != "-"]
                else:
                    level_type_abund = level_type_abund[level_type_abund[self.option("level_type")] == self.option(
                        "level_type_name")]
                    level_type_abund = level_type_abund[level_type_abund[self.option("lowest_level")] != "-"]
                level_type_abund_table = level_type_abund.groupby(self.option("lowest_level")).sum()
            else:
                level_type_abund = level_type_abund[level_type_abund[self.option("level_type")] != "-"]
                level_type_abund_table = level_type_abund.groupby(self.option("level_type")).sum()
            level_type_abund_table = level_type_abund_table.ix[list((level_type_abund_table > 0).any(axis=1))]
            if len(level_type_abund_table) < 1:  # add by zhujuan 当分类水平下数据都空时，报错并提示 20171124
                self.set_error('在所选物种/功能的分类参数下数据为空，请重新设置该参数!', code="32703901")
                raise Exception('在所选物种/功能的分类参数下数据为空，请重新设置该参数!')
            elif len(level_type_abund_table) > 300000: # add by qingchen.zhang 当分类水平下的数据超过300000时，报错并提示
                self.set_error('在所选物种/功能的分类参数下数据超过30万条，请重新设置该参数!', code="32703905")
            level_type_abund_table.to_csv(new_otu_file_path, sep="\t")

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        try:
            self.option('out_table', os.path.join(self.output_dir, "new_abund_table.xls"))
            self.logger.info(self.option('out_table').prop['path'])
        except Exception as e:
            raise Exception("输出结果文件异常——{}".format(e))

    def run(self):
        super(CreateAbundTableLefseTool, self).run()
        self.create_abund_table()
        self.set_output()
        self.end()
