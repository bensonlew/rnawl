# -*- coding: utf-8 -*-

"""群落组成分析workflow"""
import os
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
import json


class TaxonCompositionWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TaxonCompositionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "main_col", "type": "string"},
            {"name": "table", "type": "infile", "format": "meta_genomic.taxon_dir"},
            {"name": "col", "type": "string"},
            {"name": "name2id", "type": "string"},
            {"name": "level_id", "type": "int"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "others", "type": "float", "default": 0.01},
            {"name": "group_method", "type": "string", "default": ""},
            {"name": "graphic_type", "type": "string", "default": ""},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.composition_analysis = self.add_module("meta.composition.composition_analysis")
        self.filter_table = self.add_tool("metagenomic.table_select")

    def run(self):
        self.filter_table.on("end", self.run_composition_analysis)
        self.composition_analysis.on("end", self.set_db)
        self.run_filter_table()
        super(TaxonCompositionWorkflow, self).run()

    def run_composition_analysis(self):
        self.logger.info(self.filter_table.option("out_table").path)
        self.composition_analysis.set_options({
            "analysis": self.option('graphic_type'),
            "abundtable": self.filter_table.option("out_table").path,
            "group": self.option('group'),
            "add_Algorithm":  self.option('group_method'),
            "others":  self.option('others'),
        })
        self.composition_analysis.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        api_composition = self.api.api("metagenomic.composition")
        if self.option("graphic_type") in ["bubble"]:
            file_path = self.composition_analysis.output_dir.rstrip('/') + '/bubble/taxa.percents.table.xls'
            api_composition.add_composition_detail(file_path, self.option("main_id"), species_tree="",
                                                   specimen_tree="", group_method=self.option('group_method'),
                                                   level_color=self.option('level_color'), main_col="taxon_composition")
        elif self.option("graphic_type") in ["bar"]:
            file_path = self.composition_analysis.output_dir.rstrip('/') + '/bar/taxa.percents.table.xls'
            api_composition.add_composition_detail(file_path, self.option("main_id"), species_tree="",
                                                   specimen_tree="",group_method=self.option('group_method'), main_col="taxon_composition")
        elif self.option("graphic_type") in ["circos"]:
            file_path = self.composition_analysis.output_dir.rstrip('/') + '/circos/taxa.percents.table.xls'
            api_composition.add_composition_detail(file_path, self.option("main_id"), species_tree="",
                                                   specimen_tree="", group_method=self.option('group_method'), main_col="taxon_composition")
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "taxon_composition")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": "taxoncomposition",
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def run_filter_table(self):
        samples = json.loads(self.option("name2id")).keys()
        select_col = [self.option("col")] + samples
        opts = {
            "cols": json.dumps(select_col),
            "table": self.option("table").get_level(self.option("level_id"))
        }
        self.filter_table.set_options(opts)
        self.filter_table.run()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))

        if self.option("graphic_type") in ["bubble"]:
            if os.path.exists(self.output_dir + "/taxa.table.xls"):
                os.remove(self.output_dir + "/taxa.table.xls")
            os.link(self.composition_analysis.output_dir + "/bubble/taxa.table.xls", self.output_dir + "/taxa.table.xls")
            if os.path.exists(self.output_dir + "/taxa.percents.table.xls"):
                os.remove(self.output_dir + "/taxa.percents.table.xls")
            os.link(self.composition_analysis.output_dir + "/bubble/taxa.percents.table.xls", self.output_dir + "/taxa.percents.table.xls")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Bubble图结果目录", 0, "120281"],
                ["./taxa.table.xls", "xls", "物种/功能丰度结果表", 0, "120168"],
                ["./taxa.percents.table.xls", "xls", "物种/功能相对丰度结果表", 0, "120169"]
            ])
        elif self.option("graphic_type") in ["bar"]:
            if os.path.exists(self.output_dir + "/taxa.table.xls"):
                os.remove(self.output_dir + "/taxa.table.xls")
            os.link(self.composition_analysis.output_dir + "/bar/taxa.table.xls", self.output_dir + "/taxa.table.xls")
            if os.path.exists(self.output_dir + "/taxa.percents.table.xls"):
                os.remove(self.output_dir + "/taxa.percents.table.xls")
            os.link(self.composition_analysis.output_dir + "/bar/taxa.percents.table.xls", self.output_dir + "/taxa.percents.table.xls")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "柱形图结果目录", 0, "120098"],
                ["taxa.table.xls", "xls", "物种/功能丰度结果表", 0, "120168"],
                ["taxa.percents.table.xls", "xls", "物种/功能相对丰度结果表", 0, "120169"],
                ["bar.pdf", "pdf", "群落柱形图"]
            ])
        elif self.option("graphic_type") in ["circos"]:
            if os.path.exists(self.output_dir + "/taxa.table.xls"):
                os.remove(self.output_dir + "/taxa.table.xls")
            if os.path.exists(self.output_dir + "/taxa.percents.table.xls"):
                os.remove(self.output_dir + "/taxa.percents.table.xls")
            os.link(self.composition_analysis.output_dir + "/circos/taxa.table.xls", self.output_dir + "/taxa.table.xls")
            os.link(self.composition_analysis.output_dir + "/circos/taxa.percents.table.xls", self.output_dir + "/taxa.percents.table.xls")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Circos样本与物种或功能关系结果目录", 0, "120100"],
                ["taxa.table.xls", "xls", "物种/功能丰度结果表", 0, "120168"],
                ["taxa.percents.table.xls", "xls", "物种/功能相对丰度结果表", 0, "120169"]
            ])
        super(TaxonCompositionWorkflow, self).end()
