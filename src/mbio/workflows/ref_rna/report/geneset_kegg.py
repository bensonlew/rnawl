# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re


class GenesetKeggWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetKeggWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_kegg", "type": "string"},
            {"name": "kegg_table", "type": "infile", "format": "annotation.kegg.kegg_table"},
            {"name": "geneset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "add_info", "type": "string"}  # 底图颜色信息

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.kegg_class = self.add_tool("denovo_rna.express.kegg_class")

    def run(self):
        # super(GenesetKeggWorkflow, self).run()
        self.kegg_class.on("end", self.set_db)
        self.run_kegg_class()
        super(GenesetKeggWorkflow, self).run()
        # self.set_db()
        # self.end()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.ref_rna_geneset
        self.logger.info("开始进行kegg_class的导表")
        output_file = self.kegg_class.output_dir + '/kegg_stat.xls'
        pathway_file = self.kegg_class.output_dir + '/pathways'
        api_geneset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
        # api_geneset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        record_id = self.option("main_table_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = api_geneset.db["sg_geneset_kegg_class"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong")

        graph_dir = os.path.join(self.workflow_output, 'pathways')
        conn.update({"_id": record_id}, {"$set": {'graph_dir': graph_dir}}, upsert=True)
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.kegg_class.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基因集KEGG功能分类结果目录"],
        ])
        super(GenesetKeggWorkflow, self).end()

    def run_kegg_class(self):
        opts = {
            "geneset_kegg": self.option("geneset_kegg"),
            "kegg_table": self.option("kegg_table"),
            "geneset_id": self.option("geneset_id"),
            "background_links": self.option("add_info")
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()
