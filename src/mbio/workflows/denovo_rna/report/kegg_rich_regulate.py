# -*- coding: utf-8 -*-
# __author__: zengjing

"""无参转录组kegg富集、调控分析"""

from biocluster.workflow import Workflow
from biocluster.config import Config
import os
import re
import shutil
from mbio.packages.denovo_rna.express.pvalue_sort import *


class KeggRichRegulateWorkflow(Workflow):
    """
    报告中调用kegg富集分析或者调控分析时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(KeggRichRegulateWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "kegg_table", "type": "string", "default": "none"},  # 只含有基因的kegg table结果文件
            {"name": "all_list", "type": "string", "default": "none"},  # gene名字文件
            {"name": "diff_stat", "type": "string", "default": "none"},
            {"name": "correct", "type": "string", "default": "BH"},  # 多重检验校正方法
            {"name": "regulate", "type": "string", "default": "all"},
            {"name": "analysis_type", "type": "string"},
            {"name": "kegg_regulate_id", "type": "string"},
            {"name": "kegg_enrich_id", "type": "string"},
            {"name": "express_diff_id", "type": "string"},
            {"name": "sort_id", "type": "string"},
            {"name": "name", "type": "string"},
            {"name": "compare_name", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.kegg_rich = self.add_tool("denovo_rna.express.kegg_rich")
        self.kegg_regulate = self.add_tool("denovo_rna.express.kegg_regulate")
        self.rich_out = self.kegg_rich.output_dir
        self.regulate_out = self.kegg_regulate.output_dir

    def run_kegg_rich(self):
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
            "kegg_table": self.option("kegg_table"),
            "all_list": self.option("all_list"),
            "diff_list": diff_list,
            "correct": self.option("correct")
        }
        self.kegg_rich.set_options(options)
        self.kegg_rich.run()

    def run_kegg_regulate(self):
        options = {
            "kegg_table": self.option("kegg_table"),
            "diff_stat": self.option("diff_stat")
        }
        self.kegg_regulate.set_options(options)
        self.kegg_regulate.run()

    def set_output(self):
        self.logger.info("set_output")
        for f in os.listdir(self.rich_out):
            fp = os.path.join(self.rich_out, f)
            fq = os.path.join(self.output_dir, f)
            if os.path.exists(fq):
                os.remove(fq)
            shutil.copy(fp, self.output_dir)
        for f in os.listdir(self.regulate_out):
            fp = os.path.join(self.regulate_out, f)
            fq = os.path.join(self.output_dir, f)
            if os.path.exists(fq):
                try:
                    os.remove(fq)
                except:
                    os.rmdir(fq)
            try:
                shutil.copy(fp, fq)
            except:
                shutil.copytree(fp, fq)

    def pval_sort(self):
        enrich_path, regulate_path = '', ''
        for f in os.listdir(self.output_dir):
            if re.search(r".*enrichment.*$", f):
                enrich_path = os.path.join(self.output_dir, f)
            if re.search(r".*regulate.*$", f):
                regulate_path = os.path.join(self.output_dir, f)
        stat_path = os.path.join(self.work_dir, "enrich_regulate_stat.xls")
        pvalue_path = os.path.join(self.output_dir, "kegg_pvalue_sort.xls")
        if os.path.exists(enrich_path) and os.path.exists(regulate_path):
            kegg_sort_pval(enrich_path, regulate_path, stat_path, pvalue_path)
        else:
            raise Exception("output里没有kegg富集或调控的结果文件!")

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info("set_db")
        api_rich = self.api.denovo_kegg_rich
        api_regulate = self.api.denovo_kegg_regulate
        api_pval = self.api.denovo_kegg_pval_sort
        enrich_path, regulate_path, pval_path, pathway_path = '', '', '', ''
        if self.option("analysis_type") == "enrich":
            for f in os.listdir(self.rich_out):
                if re.search(r".*enrichment.*$", f):
                    enrich_path = os.path.join(self.rich_out, f)
            api_rich.add_kegg_enrich_detail(enrich_id=self.option("kegg_enrich_id"), kegg_enrich_table=enrich_path)
        if self.option("analysis_type") == "regulate":
            for f in os.listdir(self.regulate_out):
                if re.search(r".*regulate.*$", f):
                    regulate_path = os.path.join(self.regulate_out, f)
                if re.search(r".*pathways.*$", f):
                    pathway_path = os.path.join(self.regulate_out, f)
            api_regulate.add_kegg_regulate_detail(regulate_id=self.option("kegg_regulate_id"), kegg_regulate_table=regulate_path)
            api_regulate.add_kegg_regulate_pathway(regulate_id=self.option("kegg_regulate_id"), pathway_dir=pathway_path)
        if self.option("analysis_type") == "both":
            for f in os.listdir(self.output_dir):
                if re.search(r".*pvalue_sort.*$", f):
                    pval_path = os.path.join(self.output_dir, f)
            api_pval.add_pvalue_detail(sort_id=self.option("sort_id"), pval_sort=pval_path)

    def rely_fun(self):
        self.set_output()
        self.pval_sort()
        self.set_db()
        self.end()

    def run(self):
        if self.option("analysis_type") == "enrich":
            self.kegg_rich.on("end", self.set_output)
            self.kegg_rich.on("end", self.set_db)
            self.run_kegg_rich()
        if self.option("analysis_type") == "regulate":
            self.kegg_regulate.on("end", self.set_output)
            self.kegg_regulate.on("end", self.set_db)
            self.run_kegg_regulate()
        if self.option("analysis_type") == "both":
            self.on_rely([self.kegg_rich, self.kegg_regulate], self.rely_fun)
            self.run_kegg_rich()
            self.run_kegg_regulate()
        super(KeggRichRegulateWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            [r"kegg_enrichment.xls", "xls", "kegg富集分析结果"],
            [r"pathways", "dir", "kegg调控统计pathway结果图片"],
            [r"kegg_regulate_stat.xls", "xls", "kegg调控统计表"]
        ])
        super(KeggRichRegulateWorkflow, self).end()
