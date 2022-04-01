# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd

class MetabsetIpath3Workflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetIpath3Workflow, self).__init__(wsheet_object)
        options = [
            {"name": "metabset", "type": "infile", "format": "metabolome.mul_metabset"},
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ipath = self.add_tool("metabolome.metabset.ipath3")

    def run(self):
        self.ipath.on("end", self.set_db)
        self.run_ipath()
        super(MetabsetIpath3Workflow, self).run()


    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_ipath = self.api.api('metabolome.ipath')
        self.logger.info("开始进行ipath的导表")
        api_ipath.add_ipath_detail(self.option("main_table_id"), self.ipath.work_dir, self.option('metabset').path)
        # api_ipath.add_ipath_detail(self.option("main_table_id"), self.ipath.output_dir, self.option('metabset'))
        # api_proteinset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        self.end()

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        all_files = ["Metabolism.png", "Antibiotics.png", "Microbial_metabolism.png","Secondary_metabolites.png",
                     "Metabolism.svg", "Antibiotics.svg", "Microbial_metabolism.svg","Secondary_metabolites.svg","gene_ipath_input.xls"]
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            source = os.path.join(self.ipath.work_dir, fname)
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
        super(MetabsetIpath3Workflow, self).end()

    def run_ipath(self):
        opts = {
            "metabset": self.option("metabset"),
            "anno_overview": self.option("anno_overview"),
        }
        self.ipath.set_options(opts)
        self.ipath.run()
