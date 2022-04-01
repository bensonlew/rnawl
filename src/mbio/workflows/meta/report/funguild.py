# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class FunguildWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(FunguildWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile","format": "meta.otu.group_table"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "otu_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "method","type":"string","default":''},
            {"name": "others", "type" : "float", "default":0}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.funguild = self.add_tool('meta.funguild')
        self.sort_tax_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False


    def run_tax_sort_samples(self):
        abund_table = self.option("otu_table")
        self.sort_tax_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group_table'),
            "method" : self.option("method")
        })

        self.sort_tax_samples.run()


    def run_funguild(self):

        tax_abund_table =  self.sort_tax_samples.option("out_otu_table").prop['path']
        self.funguild.set_options({
            "taxon_table": tax_abund_table,
            "others" : self.option("others")
            })
        self.funguild.on("end", self.set_db)
        self.funguild.run()


    def run(self):

        self.sort_tax_samples.on("end",self.run_funguild)
        self.run_tax_sort_samples()
        super(FunguildWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        self.api_fun = self.api.api("funguild")
        sum_path = self.funguild.output_dir+"/FUNGuild_guild.txt"
        detail_path = self.funguild.output_dir+"/Funguild.txt"

        self.api_fun.add_detail(file_path=detail_path,table_id=self.option("main_id"))
        self.api_fun.add_sum(file_path=sum_path,table_id=self.option("main_id"))
        remote_dir = self._sheet.output
        self.api_fun.add_path(main_id=self.option("main_id"), path=remote_dir)
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_funguild")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "funguild",
                "interaction": 1,
                "main_table": "sg_funguild",
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.funguild.output_dir))
        save_params(self.funguild.output_dir, self.id)
        result_dir = self.add_upload_dir(self.funguild.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "FUNGuild功能预测结果目录",0,"110247"],
            ["./Funguild.txt", "xls", "FUNGuild功能预测结果表",0,"110248"],
            ["./Funguild_guild.txt", "txt", "Guild功能分类统计表",0,"110249"],
            ["./FUNGuild功能分类统计柱形图.pdf", "pdf", "FUNGuild功能分类统计柱形图", 0, ""]
        ])
        super(FunguildWorkflow, self).end()
