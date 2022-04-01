# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd

class GenesetIpathWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetIpathWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_kegg", "type": "string"},
            {"name": "kegg_table", "type": "infile", "format": "smallrna.kegg_table"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "type", "type": "string" , "default": "origin"} # 指示用origin注释还是latest注释
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ipath = self.add_tool("smallrna.geneset.ipath")

    def run(self):
        self.ipath.on("end", self.set_db)
        self.run_ipath()
        super(GenesetIpathWorkflow, self).run()


    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('smallrna.smallrna_geneset')
        self.logger.info("开始进行kegg_class的导表")
        output_file = self.ipath.output_dir + '/gene_ipath_input.xls'
        api_geneset.add_ipath_detail(self.option("main_table_id"), output_file, self.option('geneset_kegg'))
        # api_geneset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        record_id = self.option("main_table_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = api_geneset.db["sg_geneset_ipath"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        graph_dir = self.workflow_output
        conn.update({"_id": record_id}, {"$set": {'result_dir': graph_dir}}, upsert=True)
        self.end()

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        result_dir = self.add_upload_dir(self.ipath.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基因集Ipath结果目录"],
        ])
        super(GenesetIpathWorkflow, self).end()

    def run_ipath(self):
        opts = {
            "geneset_kegg": self.option("geneset_kegg"),
            "kegg_table": self.option("kegg_table"),
        }
        self.ipath.set_options(opts)
        self.ipath.run()
