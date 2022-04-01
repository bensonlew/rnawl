# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import pandas as pd

class RelationKeggpViewWorkflow(Workflow):
    """转录组关联分析 kegg通路可视化"""
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationKeggpViewWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_expdiff", "type": "string"},
            {"name": "metab_set", "type": "string"},
            {"name": "gene_expdiff", "type": "string"},
            {"name": "gene_set", "type": "string"},
            {"name": "metab_annot", "type": "string"},
            {"name": "gene_annot", "type": "string"},
            {"name": "relation_proj_type", "type": "string"},
            {"name": "relation_task_id", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "gene_group", "type": "string"},
            {"name": "metab_group", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.keggp_view = self.add_tool("metabolome.relation.keggp_view")

    def run(self):
        self.keggp_view.on("end", self.set_db)
        self.run_keggp_view()
        super(RelationKeggpViewWorkflow, self).run()

    def run_keggp_view(self):
        opts = {
            "metab_set": self.option("metab_set"),
            "gene_set": self.option("gene_set"),
            "metab_annot": self.option("metab_annot"),
            "gene_annot": self.option("gene_annot"),
            "metab_diff": self.option("metab_expdiff"),
            "gene_diff": self.option("gene_expdiff"),
            "gene_group": self.option("gene_group"),
            "metab_group": self.option("metab_group"),
            "relation_task_id": self.option("relation_task_id"),
            "relation_proj_type": self.option("relation_proj_type"),
        }
        self.keggp_view.set_options(opts)
        self.keggp_view.run()

    def set_db(self):
        api = self.api.api('metabolome.relation_keggpview')
        level_file = self.keggp_view.output_dir + "/level.txt"
        api.add_level(self.option('main_table_id'), level_file)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.keggp_view.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "代谢集KEGG功能通路结果", 0, "150037"],
            ["level.xls", "xls", "代谢集KEGG通路分类层级表", 0, "150038"],
            ["pathway_img", "", "代谢集KEGG通路图", 0, "150041"]
        ])
        super(RelationKeggpViewWorkflow, self).end()
