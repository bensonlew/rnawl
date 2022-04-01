# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

"""Venn表计算"""

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name
from comm_table import CommTableWorkflow


class VennWorkflow(CommTableWorkflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "anno_id", "type": "string"},
            {"name": "level_id", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "params", "type": "string"},
            {"name": "anno_table", "type": "infile", "format": "meta.profile"},
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            #{"name": "level_type", "type": "string", "default": ""},
            {"name": "level_type_name", "type": "string", "default": ""},
            {"name": "lowest_level", "type": "string", "default": ""},
            {"name": "abund_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "group_detail", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "method", "type": "string"},
            {"name": "main_id", "type": "string"},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1},
            #{"name": "up_path", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.venn = self.add_tool("graph.venn_table")
        self.venn_pie = self.add_tool("meta.venn_category_abund")

    def run_sort_samples(self):
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            abund_table = self.abundance.option("out_table").prop['path']
        self.logger.info(abund_table)
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option("group"),
        })
        self.sort_samples.on("end", self.run_venn)
        self.sort_samples.run()

    def run_venn(self):
        abund_table = self.sort_samples.option("out_otu_table").prop['path']
        num_lines = open(abund_table, 'r').readlines()
        if len(num_lines) < 3:  # 1个物种/功能不能计算venn
            raise OptionError("该分类下只有1个物种/功能，不能计算venn，请更换分类水平或者选择其他样品！", code="12802801")
        options = {
            "otu_table": self.sort_samples.option("out_otu_table"),
            "group_table": self.option("group"),
            "analysis_model": "mg"
        }
        self.venn.set_options(options)
        self.venn.on('end', self.run_venn_pie)
        self.venn.run()

    def run_venn_pie(self):
        venn_table = self.venn.output_dir + "/venn_table.xls"
        options = {
            "abund_file": self.sort_samples.option("out_otu_table"),
            "venn_table": venn_table
        }
        self.venn_pie.set_options(options)
        self.venn_pie.on('end', self.set_db)
        self.venn_pie.run()

    def set_db(self):
        self.logger.info("正在写入mongo数据库")
        api_venn = self.api.api("metagenomic.venn")
        venn_table_path = os.path.join(self.venn.output_dir, "venn_table.xls")
        venn_graph_path = os.path.join(self.venn.output_dir, "venn_graph.xls")
        venn_pie_path = os.path.join(self.venn_pie.output_dir, "venn_pie_abund.xls")
        api_venn.add_venn_detail(venn_table_path, self.option("main_id"))
        api_venn.add_venn_graph(venn_graph_path, self.option("main_id"))
        api_venn.add_venn_pie(venn_pie_path, self.option("main_id"))
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "venn")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = "venn_" + self.option("anno_type")
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

    def run(self):
        self.logger.info(">>>>>>>>>>>>>>>>>>>")
        self.logger.info(self._sheet)
        if self.option("abund_file").is_set:
            self.run_sort_samples()
        else:
            #self.run_get_abund_table()
            self.run_abundance(self.run_sort_samples)
            self.abundance.run()
        super(VennWorkflow, self).run()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.venn.output_dir))
        result_dir = self.add_upload_dir(self.venn.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Venn图结果目录", 0, "120172"],
            ["venn_table.xls", "xls", "共有和独有物种/功能表", 0, "120173"],
            ["venn_graph.xls", "xls", "各(组)样本物种/功能表", 0, "120174"],
            ["venn.pdf", "pdf", "venn图"]
        ])
        super(VennWorkflow, self).end()
