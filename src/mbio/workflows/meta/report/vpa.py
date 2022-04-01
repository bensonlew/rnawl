# -*- coding: utf-8 -*-
# __author__ = 'guanqing.zou'

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError
from mbio.packages.meta.common_function import envname_restore
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class VpaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VpaWorkflow, self).__init__(wsheet_object)
        options = [
             {"name": "otu_id", "type": "string"},
            {"name": "otu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            #{"name": "group_table", "type": "infile","format": "meta.otu.group_table"},
            {"name":"group","type":"infile","format": "meta.otu.group_table"},
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name":"update_info","type":"string"},
            {"name":"main_id","type":"string"},
            #{"name": "envtable", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "env_file", "type": "infile", "format": "meta.otu.otu_table"} , #guanqing.zou 20180926
            {"name": "env_group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "env_detail", "type": "string"},
            {"name": "env_labs", "type": "string"},
            {"name": "env_id", "type": "string"}

            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.sort_func_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.vpa = self.add_tool('statistical.vpa')

    def run_func_sort_samples(self):
        abund_table = self.option('otu_table')
        self.sort_func_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group'),
        })
        self.sort_func_samples.on("end", self.run_vpa)
        self.sort_func_samples.run()

    def run_vpa(self):

        otu_table = self.sort_func_samples.option("out_otu_table")
        self.vpa.set_options({
            "species_table" :  otu_table,
            "env_table" : self.option("env_file"),
            "group_table": self.option('env_group')
        })
        self.vpa.on("end", self.set_db)
        self.vpa.run()


    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_func_sort_samples()
        super(VpaWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        self.api_vpa = self.api.api("vpa")

        data_path = self.vpa.output_dir + '/env.R2adj.xls'
        graph_path = self.vpa.work_dir + '/env.plot.xls'
        self.api_vpa.add_vpa_detail(data_path,self.option("main_id"))
        png_path = self._sheet.output
        self.api_vpa.add_vpa_graph(graph_path,self.option("main_id"), png_path)
        self.output_dir = self.vpa.output_dir
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_vpa")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "vpa",
                "interaction": 1,
                "main_table": "sg_vpa",
            })
            self.figsave.run()
        else:
            self.end()

    @envname_restore
    def end(self):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "VPA结果文件目录",0,"110271"],
            ["./env.R2adj.xls", "xls", "R2adj数据",0,"110272"],
            ["./VPA分析结果图.pdf", "pdf", "VPA分析结果图", 0, ""],
            #["./env.plot.xls", "xls", "画图数据", 0],
            # ["./vpa.png", "png", "vpa图",0,"120343"],
            # ["./vpa.pdf", "pdf", "vpa图",0,"120344"]
        ])
        super(VpaWorkflow, self).end()
