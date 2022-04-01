# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


from biocluster.workflow import Workflow
import os
from mbio.packages.metaasv.common_function import calculate_abundance,link_dir


class LefseWorkflow(Workflow):
    """
    metaasv lefse分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LefseWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},##asv 表
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},##group表
            {"name": "group_detail", "type": "string"},
            {"name": "second_group_detail", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "lefse_type", "type": "string"},##以与其他的项目进行区分meta_taxon
            {"name": "lda_filter", "type": "float", "default": 2.0},
            {"name": "strict", "type": "int", "default": 0},
            {"name": "group_name", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "start_level", "type": "int", "default": 3},
            {"name": "end_level", "type": "int", "default": 7},
            {"name": "is_normalized", "type": "string", "default": "true"}##是否均一化（是否用相对丰度进行计算）
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.lefse = self.add_tool("statistical.lefse")
        self.logger.info(self.option("group_name"))

    def run_lefse(self):
        """
        运行tool lefse
        :return:
        """
        options = {
            "lefse_type":self.option("lefse_type"),
            "lefse_input": self.option("otu_file"),
            "lefse_group": self.option("group_file"),
            "lda_filter": self.option("lda_filter"),
            "strict": self.option("strict"),
            "lefse_gname": self.option("group_name"),
            "start_level": self.option("start_level"),
            "end_level": self.option("end_level"),
            "percent_abund" : self.option("is_normalized")
        }
        self.lefse.set_options(options)
        self.lefse.on("end", self.set_db)
        self.lefse.run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "LEfSe差异分析结果目录", 0, ""],
            ["./lefse_LDA.xls", "xls", "LEfSe分析lda数据表", 0, ""]
        ])
        super(LefseWorkflow, self).end()

    def set_db(self):
        """
        保存两组比较分析的结果表保存到mongo数据库中
        """
        link_dir(self.lefse.output_dir, self.output_dir)
        api_lefse = self.api.api("metaasv.lefse")
        lefse_path = self.output_dir + '/lefse_LDA.xls'
        if not os.path.isfile(lefse_path):
            self.logger.error("找不到报告文件:{}".format(lefse_path))
            self.set_error("找不到报告文件")
        api_lefse.add_species_difference_lefse_detail(file_path=lefse_path, table_id=self.option("main_id"))
        self.end()

    def run(self):
        """
        运行
        """
        self.run_lefse()
        super(LefseWorkflow, self).run()
