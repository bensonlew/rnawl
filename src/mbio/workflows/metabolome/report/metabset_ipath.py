# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd

class MetabsetIpathWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabsetIpathWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "metabset", "type": "infile", "format": "metabolome.mul_metabset"},
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ipath = self.add_tool("metabolome.metabset.ipath")

    def run(self):
        self.ipath.on("end", self.set_db)
        self.run_ipath()
        super(MetabsetIpathWorkflow, self).run()


    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_ipath = self.api.api('metabolome.ipath')
        self.logger.info("开始进行ipath的导表")
        api_ipath.add_ipath_detail(self.option("main_table_id"), self.ipath.output_dir, self.option('metabset').path)
        # api_ipath.add_ipath_detail(self.option("main_table_id"), self.ipath.output_dir, self.option('metabset'))
        # api_proteinset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        self.end()

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        all_files = ['Metabolic_pathways.svg', 'Regulatory_pathways.svg', 'Biosynthesis_of_secondary_metabolities.svg',\
                     'Metabolic_pathways.png', 'Regulatory_pathways.png', 'Biosynthesis_of_secondary_metabolities.png']
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
            ["Biosynthesis_of_secondary_metabolities.svg", "svg",
             "Biosynthesis of secondary metabolities", 0, "150045"],
            ["Biosynthesis_of_secondary_metabolities.png", "png",\
            "Biosynthesis of secondary metabolities", 0, "150045"],
            ["Metabolic_pathways.svg", "svg", "Metabolic pathways", 0, "150046"],
            ["Metabolic_pathways.png", "png", "Metabolic pathways", 0, "150046"],
            ["Regulatory_pathways.svg", "svg", "Regulatory pathways", 0, "150047"],
            ["Regulatory_pathways.png", "png", "Regulatory pathways", 0, "150047"],
        ])
        super(MetabsetIpathWorkflow, self).end()

    def run_ipath(self):
        opts = {
            "metabset": self.option("metabset"),
            "anno_overview": self.option("anno_overview"),
        }
        self.ipath.set_options(opts)
        self.ipath.run()
