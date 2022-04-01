# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import pandas as pd

class RelationIpathWorkflow(Workflow):
    """
    项目关联分析- ipath
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RelationIpathWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metab_set", "type": "string"},
            {"name": "gene_set", "type": "string"},
            {"name": "metab_annot", "type": "string"},
            {"name": "gene_annot", "type": "string"},
            {"name": "relation_task_id", "type": "string"},
            {"name": "relation_proj_type", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ipath = self.add_tool("metabolome.relation.ipath3")

    def run(self):
        self.ipath.on("end", self.set_db)
        self.run_ipath()
        super(RelationIpathWorkflow, self).run()


    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_ipath = self.api.api('metabolome.ipath')
        self.logger.info("开始进行ipath的导表")
        set_to_ipath = os.path.join(self.ipath.work_dir, "set_to_ipath.txt")
        api_ipath.add_relation_ipath(self.option("main_table_id"), set_to_ipath)
        ipath_file = os.path.join(self.ipath.output_dir, "ori_ipath_input.xls")
        api_ipath.add_relation_detail(self.option("main_table_id"), ipath_file)
        self.end()

    def end(self):
        all_files = ["Metabolism.png", "Antibiotics.png", "Microbial_metabolism.png","Secondary_metabolites.png",
                     "Metabolism.svg", "Antibiotics.svg", "Microbial_metabolism.svg","Secondary_metabolites.svg","ori_ipath_input.xls"]
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            source = os.path.join(self.ipath.output_dir, fname)
            os.link(source, link)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "iPath代谢通路分析结果"],
            ["Metabolism.svg", "svg","Metabolism", 0],
            ["Metabolism.png", "png","Metabolism", 0],
            ["Antibiotics.svg", "svg", "Antibiotics", 0 ],
            ["Antibiotics.png", "png", "Antibiotics", 0 ],
            ["Microbial_metabolism.svg", "svg", "Microbial_metabolism", 0],
            ["Microbial_metabolism.png", "png", "Microbial_metabolism", 0],
            ["Secondary_metabolites.svg", "svg", "Secondary_metabolites", 0],
            ["Secondary_metabolites.png", "png", "Secondary_metabolites", 0]
        ])
        super(RelationIpathWorkflow, self).end()

    def run_ipath(self):
        opts = {
            "metab_set": self.option("metab_set"),
            "gene_set": self.option("gene_set"),
            "metab_annot": self.option("metab_annot"),
            "gene_annot": self.option("gene_annot"),
            "relation_task_id": self.option("relation_task_id"),
            "relation_proj_type": self.option("relation_proj_type"),
        }
        self.ipath.set_options(opts)
        self.ipath.run()

