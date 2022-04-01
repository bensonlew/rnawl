# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from mbio.packages.metaasv.common_function import link_dir
from biocluster.workflow import Workflow
import os


class RankAbundanceWorkflow(Workflow):
    """
    Meta ASV Rank_abundance曲线分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RankAbundanceWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", 'format': "metaasv.otu_table"},  # 输入ASV表
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},#输入Group表
            {"name": "update_info", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "step", "type": "int", "default": 1},  # 取样步长
            {"name": "method", "type": "string", "default": "relative"},  # 计算相对丰度
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.rank_abundance = self.add_tool('metaasv.rank_abundance')
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False

    def check_option(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("otu_table").is_set:
            raise self.set_error("输入文件不存在".format(self.option("otu_table").prop['path']))

    def run(self):
        """
        运行
        :return:
        """
        self.run_rank_abundance()
        super(RankAbundanceWorkflow, self).run()

    def run_rank_abundance(self):
        """
        运行 tool
        :return:
        """
        opts = {
            "otu_table": self.option("otu_table"),
            "method": self.option("method"),
            "step" :self.option("step")
        }
        self.rank_abundance.set_options(opts)
        self.rank_abundance.on("end", self.set_db)
        self.rank_abundance.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        link_dir(self.rank_abundance.output_dir, self.output_dir)
        api_rank = self.api.api("metaasv.rank_abundance")
        rank_file = os.path.join(self.output_dir, "Rank_abundance.xls")
        if self.option("main_id"):
            main_id = self.option("main_id")
        else:
            main_id = api_rank.add_rank(asv_id="asv_id", params="params", name="Rank_abundance_Origin")
        api_rank.add_rank_detail(rank_file, main_id)
        self.end()

    def sed_files(self):
        """
        结果文件上传
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Rank_abundance曲线分析结果目录", 0, ""],
            ["./Rank_abundance.xls", "", "Rank_abundance曲线表", 0, ""]
        ])

    def end(self):
        self.sed_files()
        super(RankAbundanceWorkflow, self).end()