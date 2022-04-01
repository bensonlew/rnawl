# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

from biocluster.workflow import Workflow
import glob
import os
from bson import ObjectId
import re
from biocluster.core.exceptions import OptionError
from comm_table import CommTableWorkflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
from mbio.packages.meta.common_function import envname_restore
import json

class HeatmapCorWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(HeatmapCorWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_id", "type": "string"},
            {"name": "anno_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "env_id", "type": "string"},
            {"name": "env_labs", "type": "string"},
            {"name": "env_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "abund_file", "type": "infile", 'format': "meta.otu.otu_table"},  # 当丰度文件已存在时
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "analysis_method", "type": "string", "default": "pearsonr"},
            {"name": "env_cluster", "type": "string", "default": ""},
            {"name": "species_cluster", "type": "string", "default": ""},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "top_species", "type": "int", "default": 50},  # add new option (flit top N species)
            {"name": "main_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.correlation = self.add_tool('statistical.pearsons_correlation')

    def run_sort_samples(self):
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            abund_table = self.abundance.option("out_table").prop['path']
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option('group'),
        })
        self.sort_samples.on("end", self.run_correlation)
        self.sort_samples.run()

    def run_correlation(self):
        env_cluster = self.option("env_cluster")
        species_cluster = self.option("species_cluster")
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            abund_table = self.sort_samples.option("out_otu_table")
        self.correlation.set_options({
            "otutable": abund_table,
            "envtable": self.option('env_file'),
            "method": self.option('analysis_method'),
            "env_cluster": env_cluster,
            "species_cluster": species_cluster,
            "top_species": self.option('top_species')
            })
        self.correlation.on("end", self.set_db)
        self.correlation.run()

    def run(self):
        #self.IMPORT_REPORT_DATA = True
        #self.IMPORT_REPORT_DATA_AFTER_END = False
        if self.option("abund_file").is_set:
            self.run_sort_samples()
        else:
            #self.run_get_abund_table()
            self.run_abundance(self.run_sort_samples)
            self.abundance.run()
        super(HeatmapCorWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        species_tree = ""
        env_tree = ""
        #api_correlation = self.api.heatmap_cor
        self.api_correlation = self.api.api("metagenomic.heatmap_cor")
        corr_path = glob.glob(self.correlation.output_dir+"/*correlation*")
        pvalue_path = glob.glob(self.correlation.output_dir+"/*pvalue*")
        if self.option("species_cluster") != "":
            species_tree = self.correlation.work_dir + '/final_species_tree.tre'
        if self.option("env_cluster") != "":
            env_tree = self.correlation.work_dir + '/final_env_tree.tre'
        self.api_correlation.add_heatmap_cor_detail(corr_path[0], "correlation", self.option("main_id"),species_tree="",env_tree="")
        self.api_correlation.add_heatmap_cor_detail(pvalue_path[0], "pvalue", self.option("main_id"),
                                                    species_tree=species_tree, env_tree=env_tree)
        os.link(self.correlation.output_dir + "/pearsons_correlation.xls", self.output_dir + "/correlation.xls")
        os.link(self.correlation.output_dir + "/pearsons_pvalue.xls", self.output_dir + "/pvalue.xls")
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "heatmap_cor")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            self.logger.info(self.option('params'))
            #self.logger.info(self.option('submit_location'))
            submit_loc = json.loads(self.option('params'))["submit_location"]
            self.logger.info(submit_loc)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    @envname_restore
    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "相关性Heatmap分析结果目录", 0, "120127"],
            ["./correlation.xls", "xls", "相关系数表", 0, "120141"],
            ["./pvalue.xls", "xls", "相关系数对应p值表", 0, "120142"],
            ["./CorrHeatmap.pdf", "pdf", "相关性热图"]
        ])
        super(HeatmapCorWorkflow, self).end()
