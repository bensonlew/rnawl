# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'

import os
import gevent
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
import json


class TaxonDiffWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """reads 物种注释比较分析, 目前支持pcoa"""
        self._sheet = wsheet_object
        super(TaxonDiffWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "main_col", "type": "string"},
            {"name": "table", "type": "infile", "format": "meta_genomic.taxon_dir"},
            {"name": "col", "type": "string"},
            {"name": "name2id", "type": "string"},
            {"name": "level_id", "type": "int"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_detail", "type": "string", "default": ""},
            {"name": "test_type", "type": "string", "default": 'towgroup'},
            {"name": "correction", "type": "string", "default": "none"},
            {"name": "method", "type": "string", "default": "none"},
            {"name": "ci_meth", "type": "string", "default": ""},
            {"name": "ci_level", "type": "float", "default": 0.99},
            {"name": "tail_type", "type": "string", "default": ""},
            {"name": "paired", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "post_hoc", "type": "string", "default": ""},
            {"name": "post_hoc_level", "type": "float", "default": 0.99},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.diff = self.add_tool("statistical.metastat")
        self.filter_table = self.add_tool("metagenomic.table_select")

    def run(self):
        if self.sheet.id == "wpm_249708_54142_846362":
            gevent.spawn_later(5, self.set_db)
            super(TaxonDiffWorkflow, self).run()
            return
        if self.option("test_type") == "two_group":
            self.filter_table.on("end", self.run_two_group)
        else:
            self.filter_table.on("end", self.run_multi_group)
        self.run_filter_table()
        super(TaxonDiffWorkflow, self).run()

    def run_filter_table(self):
        samples = json.loads(self.option("name2id")).keys()
        select_col = [self.option("col")] + samples
        opts = {
            "cols": json.dumps(select_col),
            "table": self.option("table").get_level(self.option("level_id"))
        }
        self.filter_table.set_options(opts)
        self.filter_table.run()

    def run_two_group(self):
        infile = self.filter_table.option("out_table").path
        if self.option("method") == "student":
            options = {
                "student_input": infile,
                "student_group": self.option("group"),
                "student_correction": self.option("correction"),
                "student_type": self.option("tail_type"),
                "test": self.option("method"),
                "student_gname": "group_name",
                "student_coverage": self.option("ci_level")
            }
        elif self.option("method") == "mann":
            options = {
                "mann_input": infile,
                # "mann_ci": self.option("ci"),
                "mann_group": self.option("group"),
                "mann_correction": self.option("correction"),
                "mann_type": self.option("tail_type"),
                "test": self.option("method"),
                "mann_gname": "group_name",
                "mann_coverage": self.option("ci_level")
            }
        elif self.option("method") == "signal":
            options = {
                "signal_input": infile,
                "signal_group": self.option("group"),
                "signal_correction": self.option("correction"),
                "signal_type": self.option("tail_type"),
                "test": self.option("method"),
                "signal_gname": "group_name",
                "signal_pair_file": self.option("paired"),
                "signal_coverage": self.option("ci_level")
            }
        else:
            options = {
                "welch_input": infile,
                "welch_group": self.option("group"),
                "welch_correction": self.option("correction"),
                "welch_type": self.option("tail_type"),
                "test": self.option("method"),
                "welch_gname": "group_name",
                "welch_coverage": self.option("ci_level")
            }
        self.diff.set_options(options)
        self.diff.on("end", self.set_db)
        self.diff.run()

    def run_multi_group(self):
        infile = self.filter_table.option("out_table").path
        if self.option("method") == "anova":
            options = {
                "anova_input": infile,
                "anova_group": self.option("group"),
                "anova_correction": self.option("correction"),
                "test": self.option("method"),
                "anova_gname": "group_name",
                "anova_methor": self.option("post_hoc"),
                "anova_coverage": self.option("post_hoc_level")
            }
        else:
            options = {
                "kru_H_input": infile,
                "kru_H_group": self.option("group"),
                "kru_H_correction": self.option("correction"),
                "test": self.option("method"),
                "kru_H_gname": "group_name",
                "kru_H_methor": self.option("post_hoc"),
                "kru_H_coverage": self.option("post_hoc_level")
            }
        self.diff.set_options(options)
        self.diff.on("end", self.set_db)
        self.diff.run()

    def set_db(self):
        self.output_dir = self.diff.output_dir
        api_diff = self.api.api("metagenomic.metastat")
        api_diff.change_main_col("taxon_diff")
        if self.option("test_type") == "two_group":
            stat_path = self.output_dir + '/' + self.option("method") + '_result.xls'
            boxfile_path = self.output_dir + '/' + self.option("method") + '_boxfile.xls'
            ci_path = self.output_dir + '/' + self.option("method") + '_CI.xls'
            bar_path = self.diff.work_dir + '/' + self.option("method") + '_plot_group_bar.xls'
            if not os.path.isfile(stat_path):
                self.set_error("找不到报告文件:%s", variables=(stat_path), code="12805201")
            if not os.path.isfile(boxfile_path):
                self.set_error("找不到报告文件:%s", variables=(boxfile_path), code="12805202")
            if not os.path.isfile(ci_path):
                self.set_error("找不到报告文件:%s", variables=(ci_path), code="12805203")
            api_diff.add_metastat_detail(statfile=stat_path, cifiles=[ci_path], check_type='two_group',
                                         metastat_id=self.option('main_id'), main_col="taxon_diff")
            api_diff.add_metastat_bar(bar_path, self.option('main_id'))
            api_diff.add_metastat_box(boxfile_path, self.option('main_id'))
            api_diff.update_metastat(self.option('main_id'), stat_path, ci_path, 'two_group')
        else:
            stat_path = self.output_dir + '/' + self.option("method") + '_result.xls'
            boxfile_path = self.output_dir + '/' + self.option("method") + '_boxfile.xls'
            bar_path = self.diff.work_dir + '/' + self.option("method") + '_plot_group_bar.xls'
            cifiles = []
            for r, d, f in os.walk(self.output_dir):
                for i in f:
                    if self.option("post_hoc") in i:
                        ci_path = r + '/' + i
                        if not os.path.isfile(ci_path):
                            self.set_error("找不到报告文件:%s", variables=(ci_path), code="12802003")
                        cifiles.append(ci_path)
            api_diff.add_metastat_detail(statfile=stat_path, cifiles=cifiles, check_type='multiple',
                                         metastat_id=self.option('main_id'), posthoc=self.option("post_hoc"),
                                         main_col="taxon_diff")
            api_diff.add_metastat_bar(bar_path, self.option('main_id'))
            api_diff.add_metastat_box(boxfile_path, self.option('main_id'))
            api_diff.update_detail(self.option('main_id'))

        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "taxon_diff")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            if self.option("test_type") == "two_group":
                submit_loc = "taxondiff_two_group"
            else:
                submit_loc = "taxondiff_multi_group"
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

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("test_type") == "two_group":
            result_dir.add_relpath_rules([
                [".", "", "两组比较结果目录", 0, "120185"],
                ["DiffTwoGroup.pdf", "pdf", "两组差异检验柱形图"]
            ])
            result_dir.add_regexp_rules([
                [r".*_result\.xls", "xls", "差异显著性比较结果表，包括均值，标准差，p值", 0, "120186"],
                [r".*_CI\.xls", "xls", "差异显著性比较,两组比较置信区间值以及效应量", 0, "120188"],
                [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值", 0, "120187"]

            ])
        else:
            result_dir.add_relpath_rules([
                [".", "", "多组比较结果目录", 0, "120181"],
                ["DiffMulGroup.pdf", "pdf", "多组差异检验柱形图"]
            ])
            result_dir.add_regexp_rules([
                [r".*_result\.xls", "xls", "差异显著性比较结果表，包括均值，标准差，p值", 0, "120183"],
                [r".*(-).*", "xls", "差异显著性比较多组比较的posthoc检验比较的结果，包含置信区间，效果量，p值", 0, "120182"],
                [r".*_boxfile\.xls", "xls", "差异显著性比较箱线图数据", 0, "120184"]
            ])
        super(TaxonDiffWorkflow, self).end()
