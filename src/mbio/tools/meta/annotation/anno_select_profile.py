# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last modify date: 20171012

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class AnnoSelectProfileAgent(Agent):
    """
    宏基因组注释交互分析基因丰度表筛选
    """

    def __init__(self, parent):
        super(AnnoSelectProfileAgent, self).__init__(parent)
        options = [
            {"name": "geneset_list", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "samples", "type": "string","default": "all"},
            {"name": "sel_gene_profile", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("geneset_list").is_set:
            raise OptionError("必须设置基因list文件")
        if not self.option("gene_profile").is_set:
            raise OptionError("必须设置基因丰度文件")
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        super(AnnoSelectProfileAgent, self).end()


class AnnoSelectProfileTool(Tool):
    def __init__(self, config):
        super(AnnoSelectProfileTool, self).__init__(config)

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoSelectProfileTool, self).run()
        self.run_select()
        self.set_output()
        self.end()

    def run_select(self):
        self.logger.info("start gene profile select")
        gene_file = self.option("geneset_list").prop["path"]
        profile = self.option("gene_profile").prop["path"]
        samples = self.option("samples")
        self.logger.info(gene_file)
        self.logger.info(profile)
        self.logger.info(samples)
        gene_list = pd.read_table(gene_file, sep='\t', header=0)
        gene_list = pd.DataFrame(gene_list)
        genes = gene_list['GeneID']
        profiletable = pd.read_table(profile, sep='\t', header=0)
        df_profile = pd.DataFrame(profiletable)
        try:
            if samples != "all":
                choose_name = "GeneID," + samples
                choose_names = choose_name.split(",")
                sample_list = samples.split(",")
                self.logger.info(sample_list)
                select = df_profile[df_profile['GeneID'].isin(genes)][choose_names]
                self.logger.info("start Total")
                select['Total']= select.loc[:,sample_list].apply(lambda x: x.sum(), axis = 1 )
            else:
                select = df_profile[df_profile['GeneID'].isin(genes)]
            self.logger.info(self.output_dir + "/select_gene_profile.xls")
            select.to_csv(self.output_dir + "/select_gene_profile.xls", sep = "\t", index=False)
        except Exception as e:
            raise Exception("创建新基因丰度表失败——{}".format(e))

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        try:
            self.option("sel_gene_profile", self.output_dir + "/select_gene_profile.xls")
            self.logger.info("设置输出结果文件成功")
        except Exception as e:
            raise Exception("输出结果文件异常——{}".format(e))
