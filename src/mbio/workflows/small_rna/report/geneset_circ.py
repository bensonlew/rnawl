# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd

class GenesetCircWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetCircWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "enrich_type", "type": "string", "default": "kegg"},
            {"name": "go_type", "type": "string", "deafult":"ALL"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "diff_id", "type": "string"},
            {"name": "p_thre", "type": "float"},
            {"name": "padj_thre", "type": "float"},
            {"name": "anno_num_thre", "type": "int"},
            {"name": "enrich_id", "type": "string"},
            {"name": "compare_group", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "diff_fc", "type": "infile", "format": "smallrna.common"},
            {"name": "enrich_table", "type": "infile", "format": "smallrna.common"},
            {"name": "gene_detail", "type": "string", "default": ""},
            {"name": "anno_list", "type": "int", "default": ""},
            {"name": "geneset_id", "type": "string"},
            {"name": "type", "type": "string"} # 指示用origin注释还是latest注释
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.circ = self.add_tool("smallrna.geneset.enrich2circ")

    def run(self):
        self.circ.on("end", self.set_db)
        self.run_circ()
        super(GenesetCircWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('smallrna.smallrna_geneset')
        self.logger.info("开始进行富集弦图的导表")
        if self.option("enrich_type") == "GO":
            enrich_circ_file = self.circ.output_dir + '/go_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/go_enrich_detail.table'
        elif self.option("enrich_type") == "KEGG":
            enrich_circ_file = self.circ.output_dir + '/kegg_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/kegg_enrich_detail.table'
        enrich_zscore_file = self.circ.output_dir + '/enrich_zscore'
        api_geneset.add_circ_graph(self.option("main_table_id"), enrich_circ_file, self.option('enrich_type'))
        api_geneset.add_circ_detail(self.option("main_table_id"), enrich_detail_file, self.option('enrich_type'))
        api_geneset.update_circ_main(self.option("main_table_id"), enrich_zscore_file, self.option('enrich_type'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.circ.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基因集Circ结果目录"],
        ])
        super(GenesetCircWorkflow, self).end()

    def run_circ(self):
        opts = {
            "enrich_type" : self.option("enrich_type"),
            "enrich_table" : self.option("enrich_table"),
            "fc_table" : self.option("diff_fc"),
            "p_thre" : self.option("p_thre"),
            "padj_thre" : self.option("padj_thre"),
            "anno_num_thre" : self.option("anno_num_thre"),
            "anno_list" : self.option("anno_list"),
            "gene_detail" : self.option("gene_detail")
        }
        self.circ.set_options(opts)
        self.circ.run()
