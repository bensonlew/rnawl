# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class RecalAnnoAbuModule(Module):
    """
    宏基因组基础注释和个性化注释交互分析使用，筛选后gene_anno表重新计算各层级丰度等
    """
    def __init__(self, work_id):
        super(RecalAnnoAbuModule, self).__init__(work_id)
        options = [
            {"name": "gene_anno_table", "type": "infile", "format": "sequence.profile_table", "required":"True"},
            {"name": "gene_profile", "type": "infile", "format": "sequence.profile_table", "required":"True"},
            {"name": "go1234level_out", "type": "infile", "format": "sequence.profile_table"},
            {"name": "database", "type": "string"},  # nr、cog、kegg、cazy、ardb、card、vfdb
            {"name": "origin_anno", "type": "infile", "format": "sequence.profile_table"},  # kegg时使用
            {"name": "xml_file", "type": "infile", "format": "sequence.profile_table"},  # kegg时使用
            {"name": "anno_result_dir", "type": "outfile", "format": "annotation.mg_anno_dir"},
            {"name": "vfdb_type", "type": "string", "default": "core"},  # vfdb时使用
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "result_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            {"name": "vfdb_type", "type": "string", "default": "all"},  # vfdb使用
            {"name": "task_id", "type": "string", "default": ""}, # 用于区分新老任务
        ]
        self.add_option(options)
        self.basic_database_tool = self.add_tool("meta.annotation.recal_anno_abu")  # nr、cog、kegg、cazy、ardb、card、vfdb

    def check_options(self):
        if not self.option("database"):
            raise OptionError("必须设置数据库名称database!", code="21202701")
        return True

    def set_output(self, event):
        self.end()

    def run(self):
        super(RecalAnnoAbuModule, self).run()
        self.basic_anno = ["nr", "cog", "kegg", "cazy", "ardb", "card", "vfdb"]
        self.personal_anno = ["go", "phi", "mvirdb", "tcdb", "qs", "pfam", "sec", "t3ss", "probio", "p450"]
        if self.option("database") in self.basic_anno:
            self.run_basic_anno()
        elif self.option("database") in self.personal_anno:
            self.run_personal_anno()
        else:
            self.set_error("database name %s is wrong", variables=(self.option("database")), code="21202701")

    def run_basic_anno(self):
        self.logger.info("basic_anno recalculate!")
        options = {
            "gene_anno_table": self.option("gene_anno_table"),
            "gene_profile": self.option("gene_profile"),
            "database": self.option("database"),
            "task_id": self.option("task_id")
        }
        if self.option("database") == "kegg":
            # options["xml_file"] = self.option("xml_file")
            options["origin_anno"] = self.option("origin_anno")
        if self.option("database") == "vfdb":
            options["vfdb_type"] = self.option("vfdb_type")
        if self.option("group").is_set:
            options["group"] = self.option("group")
        self.basic_database_tool.set_options(options)
        self.basic_database_tool.on("end", self.end)
        self.basic_database_tool.run()

    def run_personal_anno(self):
        self.logger.info("personal_anno recalculate!")
        gene_anno_table = self.option("gene_anno_table").prop["path"]
        gene_profile = self.option("gene_profile").prop["path"]
        options = {}
        if self.option("database") == "go":
            self.personal_anno_tool = self.add_module("annotation.anno_go_stat")
            options["go1234level_out"] = self.option("go1234level_out").prop['path']
            options["reads_profile_table"] = gene_profile
            options["gene_anno"] = gene_anno_table
        if self.option("database") == "qs":
            self.personal_anno_tool = self.add_tool("annotation.qs_anno_stat")
            if self.option("group").is_set:
                options["group_table"] = self.option("group")
            options["qs_anno_result"] = gene_anno_table
        if self.option("database") == "probio":
            self.personal_anno_tool = self.add_tool("annotation.probio_abun")
            options["probio_anno"] = gene_anno_table
        if self.option("database") == "pfam":
            self.personal_anno_tool = self.add_tool("annotation.pfam_anno_stat")
            options["pfam_anno_table"] = gene_anno_table
        if self.option("database") == "p450":
            self.personal_anno_tool = self.add_tool("annotation.cyps_anno_stat")
            options["cyps_anno_table"] = gene_anno_table
        if self.option("database") == "tcdb":
            self.personal_anno_tool = self.add_tool("annotation.tcdb_anno_stat")
            options["tcdb_anno_table"] = gene_anno_table
        if self.option("database") == "mvirdb":
            self.personal_anno_tool = self.add_tool("annotation.mvirdb_anno_stat")
            options["mvirdb_anno_table"] = gene_anno_table
        if self.option("database") == "phi":
            self.personal_anno_tool = self.add_tool("annotation.phi_anno_stat")
            options["phi_anno_table"] = gene_anno_table
        options["reads_profile_table"] = gene_profile
        self.personal_anno_tool.set_options(options)
        self.personal_anno_tool.on("end", self.end)
        self.personal_anno_tool.run()

    def end(self):
        if self.option("database") in self.basic_anno:
            self.output_dir = self.basic_database_tool.output_dir
        elif self.option("database") in self.personal_anno:
            self.output_dir = self.personal_anno_tool.output_dir
        self.option("result_dir",self.output_dir)
        super(RecalAnnoAbuModule, self).end()
