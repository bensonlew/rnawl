# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

"""无参转录组go富集分析/go调控分析"""

from biocluster.workflow import Workflow
from biocluster.config import Config
import os
import re
import shutil
from mbio.packages.denovo_rna.express.pvalue_sort import *


class GoEnrichRegulateWorkflow(Workflow):
    """
    交互分析时调用go富集分析/go调控分析时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GoEnrichRegulateWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "all_list", "type": "string", "default": "none"},
            {"name": "go_list", "type": "string", "default": "none"},
            {"name": "diff_stat", "type": "string", "default": "none"},
            {"name": "go2level", "type": "string", "default": "none"},
            {"name": "pval", "type": "string", "default": "0.05"},                        # 显著性水平
            {"name": "method", "type": "string", "default": "bonferroni,sidak,holm,fdr"}, # 多重校正方法
            {"name": "regulate", "type": "string", "default": "all"},
            {"name": "go_enrich_id", "type": "string"},
            {"name": "go_regulate_id", "type": "string"},
            {"name": "sort_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "analysis_type", "type": "string"},    # 分析类型(enrich/regulate/stat)
            {"name": "express_diff_id", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "compare_name", "type": "string"},
            {"name": "submit_location", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.go_enrich = self.add_tool("denovo_rna.express.go_enrich")
        self.go_regulate = self.add_tool("denovo_rna.express.go_regulate")
        self.enrich_output = self.go_enrich.output_dir
        self.regulate_output = self.go_regulate.output_dir

    def run_go_enrich(self):
        diff_list = os.path.join(self.work_dir, "{}_vs_{}.list".format(self.option("name"), self.option("compare_name")))
        with open(self.option("diff_stat"), "rb") as f, open(diff_list, "wb") as w:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                if self.option("regulate") == "all":
                    if float(line[6]) < 0.05:
                        w.write(line[0] + '\n')
                if self.option("regulate") == "up":
                    if float(line[6]) < 0.05 and line[9] == "up":
                        w.write(line[0] + '\n')
                if self.option("regulate") == "down":
                    if float(line[6]) < 0.05 and line[9] == "down":
                        w.write(line[0] + '\n')
        options = {
            "diff_list": diff_list,
            "all_list": self.option("all_list"),
            "go_list": self.option("go_list"),
            "pval": self.option("pval"),
            "method": self.option("method"),
        }
        self.go_enrich.set_options(options)
        self.go_enrich.run()

    def run_go_regulate(self):
        options = {
            "diff_stat": self.option("diff_stat"),
            "go_level_2": self.option("go2level"),
        }
        self.go_regulate.set_options(options)
        self.go_regulate.run()

    def set_output(self):
        self.logger.info("set_output")
        for f_name in os.listdir(self.enrich_output):
            fp = os.path.join(self.enrich_output, f_name)
            fq = os.path.join(self.output_dir, f_name)
            if os.path.exists(fq):
                os.remove(fq)
            shutil.copy(fp, self.output_dir)
        for f_name in os.listdir(self.regulate_output):
            fp = os.path.join(self.regulate_output, f_name)
            fq = os.path.join(self.output_dir, f_name)
            if os.path.exists(fq):
                os.remove(fq)
            shutil.copy(fp, self.output_dir)

    def pval_sort(self):
        method = self.option("method")
        enrich_path, regulate_path = '', ''
        for f in os.listdir(self.output_dir):
            if re.search(r"go_enrich.*$", f):
                enrich_path = os.path.join(self.output_dir, f)
            if re.search(r"GO_regulate.*$", f):
                regulate_path = os.path.join(self.output_dir, f)
        stat_path = os.path.join(self.output_dir, 'enrich_regulate_stat.xls')
        pvalue_path = os.path.join(self.output_dir, 'go_pvalue_sort.xls')
        if os.path.exists(enrich_path) and os.path.exists(regulate_path):
            go_sort_pval(enrich_path=enrich_path, regulate_path=regulate_path, stat_path=stat_path, pvalue_path=pvalue_path)

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("set_db")
        api_enrich = self.api.denovo_go_enrich
        api_regulate = self.api.denovo_go_regulate
        api_pvalue = self.api.denovo_go_pval_sort
        go_regulate_dir, go_enrich_dir, go_sort_dir, go_graph_dir = '', '', '', ''
        if self.option("analysis_type") == "enrich":
            for f in os.listdir(self.enrich_output):
                if re.search(r"go_enrich_.*$", f):
                    go_enrich_dir = self.enrich_output + '/' + f
                if re.search(r"go_lineage.*$", f):
                    go_graph_dir = self.enrich_output + '/' + f
            api_enrich.add_go_enrich_detail(go_enrich_id=self.option("go_enrich_id"), go_enrich_dir=go_enrich_dir)
            if os.path.exists(go_graph_dir):
                api_enrich.update_directed_graph(go_enrich_id=self.option("go_enrich_id"), go_graph_dir=go_graph_dir)
        if self.option("analysis_type") == "regulate":
            for f in os.listdir(self.regulate_output):
                if re.search(r"GO_regulate.xls$", f):
                    go_regulate_dir = self.regulate_output + '/' + f
            api_regulate.add_go_regulate_detail(go_regulate_id=self.option("go_regulate_id"), go_regulate_dir=go_regulate_dir)
        if self.option("analysis_type") == "both":
            for f in os.listdir(self.output_dir):
                if re.search(r"go_pvalue_sort.xls$", f):
                    go_sort_dir = self.output_dir + '/' + f
            api_pvalue.add_pvalue_detail(sort_id=self.option("sort_id"), method=self.option("method"), pval_sort=go_sort_dir)

    def rely_fun(self):
        self.set_output()
        self.pval_sort()
        self.set_db()
        self.end()

    def run(self):
        if self.option("analysis_type") == "enrich":
            self.go_enrich.on("end", self.set_output)
            self.go_enrich.on("end", self.set_db)
            self.run_go_enrich()
        if self.option("analysis_type") == "regulate":
            self.go_regulate.on("end", self.set_output)
            self.go_regulate.on("end", self.set_db)
            self.run_go_regulate()
        if self.option("analysis_type") == "both":
            self.on_rely([self.go_enrich, self.go_regulate], self.rely_fun)
            self.run_go_enrich()
            self.run_go_regulate()
        super(GoEnrichRegulateWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        relpath = [
            [".", "", "结果输出目录"],
            [r"go_enrich_.*", "xls", "go富集结果文件"],
            [r"go_lineage.*", "png", "go富集有向无环图"],
            ["GO_regulate", "xls", "基因上下调在GO2level层级分布情况表"],
            ["go_pvalue_sort", "xls", "go富集调控统计排序表"]
        ]
        result_dir.add_regexp_rules([
            [r"go_enrich_.*", "xls", "go富集结果文件"],
            [r"GO_regulate.xls", "xls", "上下调基因go注释"],
            [r"go_lineage.*", "png", "go富集有向无环图"],
            ["go_pvalue_sort", "xls", "go富集调控统计排序表"]
        ])
        result_dir.add_relpath_rules(relpath)
        super(GoEnrichRegulateWorkflow, self).end()
