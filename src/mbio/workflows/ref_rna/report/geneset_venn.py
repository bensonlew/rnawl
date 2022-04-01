# -*- coding: utf-8 -*-
# __author__ = 'konghualei, 20170421'
import web
import json
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import pandas as pd
import shutil
import re,os
from biocluster.workflow import Workflow
from mbio.packages.ref_rna.express.genesetVenn import ExpressVenn

class GenesetVennWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetVennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"geneset_file","type":"string"},
            {"name":"type","type":"string"}, #基因还是转录本
            {"name":"update_info","type":"string"},
            {"name":"geneset_venn_id","type":"string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.geneset_venn = self.add_tool("rna.geneset_venn")
        
    def run_geneset_venn(self):
        """
        基因集venn图分析
        """
        opts = {
            "fpkm":self.option("geneset_file")
        }
        self.geneset_venn.set_options(opts)
        self.geneset_venn.on("end", self.set_db)
        self.geneset_venn.run()
        
    def set_db(self):
        for files in os.listdir(self.geneset_venn.output_dir):
            shutil.copy2(self.geneset_venn.output_dir+"/"+files,self.output_dir+"/"+files)
        api_geneset_venn = self.api.refrna_corr_express
        venn_table = self.geneset_venn.output_dir + "/venn_table.xls"
        venn_id = self.option("geneset_venn_id")
        api_geneset_venn.add_venn_detail(venn_table, venn_id, project = 'ref',analysis_name="geneset")
        self.logger.info("venn_table上传完毕！")
        venn_graph = self.geneset_venn.output_dir + "/venn_graph.xls"
        api_geneset_venn.add_venn_graph(venn_graph, venn_id, project='ref',analysis_name="geneset")
        self.logger.info("venn_graph上传完毕！")
        self.end()
    
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
                ["./venn_graph.xls","xls","venn_graph结果目录"],
                ["./venn_table.xls","xls","venn_table结果目录"]
        ])
        super(GenesetVennWorkflow, self).end()
    
    def run(self):
        self.run_geneset_venn()
        super(GenesetVennWorkflow, self).run()

    def end(self):
        output1_dir = self.geneset_venn.output_dir
        result = self.add_upload_dir(output1_dir)
        result.add_relpath_rules([[".", "", "基因集Venn分析结果文件"], ])
        super(GenesetVennWorkflow, self).end()
