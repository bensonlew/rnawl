# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2017.12.27

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
import pandas as pd
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class ContributeWorkflow(CommTableWorkflow):
    """
    宏基因组Contribute注释交互分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ContributeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "gene_list", "type": "infile", "format": "sequence.profile_table"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "fun_anno_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "nr_anno_file", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "samples", "type": "string"},
            {"name": "nr_level", "type": "string"},
            {"name": "fun_level", "type": "string"},
            {"name": "top_tax", "type": "int", "default": 10},
            {"name": "top_fun", "type": "int", "default": 10},
            {"name": "update_info", "type": "string"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "group_all", "type": "string", "default": "F"},
            {"name": "group_method", "type": "int","default": 1},  # 分组求和方法：无：0，求和：1，均值：2，中位数：3
            #{"name": "params", "type": "string"},
            #{"name": "name", "type": "string"},
            #{"name": "group_detail", "type": "string", "default": ""},
            {"name": "main_table_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.gene_profile_select = self.add_tool("sequence.select_table")  # 基因profile筛选
        self.contribute_tool = self.add_tool("meta.association_model.contribute")  # 贡献度分析各样本计算
        self.run_group_tools = []
        self.tax_fun_tool = ""
        self.fun_tax_tool = ""

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start contribute!")
        if self.option("method") == "ppm" and self.option("clean_stat").is_set:
            self.run_cal_ppm(rely=self.run_profile_select, is_run=True)
        else:
            self.run_profile_select()
        super(ContributeWorkflow, self).run()

    def run_profile_select(self):
        self.logger.info("start profile!")
        gene_list = self.option("gene_list")
        if self.option("method") == "ppm" and self.option("clean_stat").is_set:
            gene_profile = self.cal_ppm.option("out_table").prop["path"]
        else:
            gene_profile = self.option("gene_profile")
        self.logger.info(gene_list)
        self.logger.info(gene_profile)
        options = {
            "select_genes": gene_list,
            "origin_table": gene_profile,
            "samples": self.option("samples")
        }
        self.gene_profile_select.set_options(options)
        self.gene_profile_select.on('end', self.run_contribute)
        self.gene_profile_select.run()

    def run_contribute(self):
        self.logger.info(self.option("top_tax"))
        sel_gene_profile = self.gene_profile_select.option("select_table")
        fun_anno_file = self.option("fun_anno_file")
        options = {
            "taxon_file": self.option("nr_anno_file"),
            "function_file": fun_anno_file,
            "gene_profile": self.gene_profile_select.option("select_table"),
            "tax_level": self.option("nr_level"),
            "fun_level": self.option("fun_level"),
            "top_tax": self.option("top_tax"),
            "top_fun": self.option("top_fun")
        }
        self.contribute_tool.set_options(options)
        # self.select_tools.append(anno_select)
        self.contribute_tool.on('end', self.run_group)
        self.contribute_tool.run()

    def run_group(self):
        files = os.listdir(self.contribute_tool.output_dir)
        self.logger.info(self.option("group_table"))
        gene_list = self.option("gene_list")
        for file in files:
            self.logger.info(self.option("group_table"))
            file_path = os.path.join(self.contribute_tool.output_dir, file)
            with open(file_path, 'r') as r:
                header = r.readline().strip().split('\t')
            if file == "Taxon_function_abundance.xls":
                self.tax_fun_tool = self.add_tool("sequence.select_table")
                options = {
                    "origin_table": file_path,
                    "group": self.option("group_table"),
                    "group_all": self.option("group_all"),
                    "group_method" : self.option("group_method"),
                    "select_columns": "{},{}".format(header[0], header[1])
                }
                self.tax_fun_tool.set_options(options)
                self.run_group_tools.append(self.tax_fun_tool)
            elif file == "Function_taxon_abundance.xls":
                self.fun_tax_tool = self.add_tool("sequence.select_table")
                options = {
                    "origin_table": file_path,
                    "group": self.option("group_table"),
                    "group_all": self.option("group_all"),
                    "group_method" : self.option("group_method"),
                    "select_columns": "{},{}".format(header[0], header[1])
                }
                self.fun_tax_tool.set_options(options)
                self.run_group_tools.append(self.fun_tax_tool)
        self.on_rely(self.run_group_tools, self.set_db)
        for tool in self.run_group_tools:
            tool.run()

    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        api_contribute = self.api.api("metagenomic.contribute")
        #specimen = self.option("samples")
        tax_fun_old = os.path.join(self.tax_fun_tool.output_dir , "select_table.xls")
        fun_tax_old = os.path.join(self.fun_tax_tool.output_dir , "select_table.xls")
        tax_fun_file = os.path.join(self.output_dir , "Taxon_function_abundance.xls")
        fun_tax_file = os.path.join(self.output_dir , "Function_taxon_abundance.xls")
        if os.path.exists(tax_fun_file):
            os.remove(tax_fun_file)
        os.link(tax_fun_old, tax_fun_file)
        if os.path.exists(fun_tax_file):
            os.remove(fun_tax_file)
        os.link(fun_tax_old, fun_tax_file)
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_table_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12801301")
        self.logger.info(main_id)
        #api_contribute.add_contribute(main_id, contribute_dir, "core", update_main=False)
        api_contribute.add_contribute_detail(main_id, tax_fun_file, fun_tax_file)
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_table_id"), "contribute")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_table_id"), "contribute")
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_table_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "贡献度分析结果目录", 0, "120211"],
            ["Function_taxon_abundance.xls", "xls", "功能的物种贡献度结果表", 0, "120212"],
            ["Taxon_function_abundance.xls", "xls", "物种的功能贡献度结果度", 0, "120213"],
            ["Contribute_species_heatmap.pdf", "pdf", "物种贡献度heatmap图"],
            ["Contribute_function_heatmap.pdf", "pdf", "功能贡献度heatmap图"],
            ["Contribute_species_bar.pdf", "pdf", "物种贡献度bar图"],
            ["Contribute_function_bar.pdf", "pdf", "功能贡献度bar图"],
        ])
        super(ContributeWorkflow, self).end()
