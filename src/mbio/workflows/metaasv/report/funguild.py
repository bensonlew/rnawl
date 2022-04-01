# -*- coding: utf-8 -*-
##__author__: qingchen.zhang
from biocluster.workflow import Workflow
from mbio.packages.metaasv.common_function import link_dir


class FunguildWorkflow(Workflow):
    """
    metaasv FunGuild预测分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FunguildWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},##OTU 表
            {"name": "level", "type": "int", "default": 9},##分类学水平
            {"name": "group_table", "type": "infile","format": "meta.otu.group_table"},## groupfile
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "method","type":"string","default":''},##分组方法计算方法
            {"name": "others", "type" : "float", "default":0}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.funguild = self.add_tool('metaasv.funguild')
        self.sort_tax_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False

    def run_tax_sort_samples(self):
        """
        排序和按分组对样本进行合并
        """
        abund_table = self.option("otu_table")
        self.sort_tax_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group_table'),
            "method" : self.option("method")
        })
        self.sort_tax_samples.run()

    def run_funguild(self):
        """
        funguild 预测分析
        """
        tax_abund_table =  self.sort_tax_samples.option("out_otu_table").prop['path']
        self.funguild.set_options({
            "taxon_table": tax_abund_table,
            "others" : self.option("others")
            })
        self.funguild.on("end", self.set_db)
        self.funguild.run()


    def run(self):
        """
        运行
        """
        self.sort_tax_samples.on("end",self.run_funguild)
        self.run_tax_sort_samples()
        super(FunguildWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中和连接结果文件
        """
        self.logger.info("正在写入mongo数据库")
        self.api_fun = self.api.api("metaasv.funguild")
        link_dir(self.funguild.output_dir, self.output_dir)
        sum_path = self.output_dir+"/FUNGuild_guild.txt"
        detail_path = self.output_dir+"/Funguild.txt"
        self.api_fun.add_detail(file_path=detail_path,table_id=self.option("main_id"))
        self.api_fun.add_sum(file_path=sum_path,table_id=self.option("main_id"))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.funguild.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "FUNGuild功能预测结果目录",0,""],
            ["./Funguild.txt", "xls", "FUNGuild功能预测结果表",0,""],
            ["./Funguild_guild.txt", "txt", "Guild功能分类统计表",0,""]
        ])
        super(FunguildWorkflow, self).end()
