# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import tarfile
import shutil
import pandas as pd
from mbio.packages.prok_rna.chart import Chart
import glob


class GenesetKeggWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetKeggWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_kegg", "type": "infile", "format": "prok_rna.common"},
            {"name": "kegg_table", "type": "infile", "format": "prok_rna.kegg_table"},
            {"name": "kegg_table_2", "type": "infile", "format": "denovo_rna_v2.vcf"},
            {"name": "geneset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "add_info", "type": "string"},  # 底图颜色信息
            {"name": "type", "type": "string"}, # 指示用origin注释还是latest注释
            {'name': 'kegg_version', 'type': 'string', 'default': "2017"},
            {"name": "task_id", "type": "string"},
            {"name": "source", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.kegg_class = self.add_tool("prok_rna.geneset.kegg_class")

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
        api_geneset = self.api.api('prok_rna.geneset')

        kegg_stat = self.kegg_class.output_dir + '/kegg_stat.xls'
        pathway_class = self.option("kegg_table_2").prop['path']
        pathway_file = self.kegg_class.output_dir + '/pathways'

        self.logger.info("开始进行kegg_class的导表")
        # api_geneset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
        #api_geneset.add_kegg_regulate_new(self.option("main_table_id"), self.option("geneset_id"), output_file, self.option("kegg_table_2").prop['path'], self.work_dir)
        api_geneset.add_kegg_regulate_new2(
            self.option("main_table_id"),
            self.option("geneset_kegg").prop['path'],
            kegg_stat,
            self.option("kegg_table_2").prop['path'],
            source=self.option('source')
        )
        api_geneset.add_kegg_regulate_pic(
            self.option("main_table_id"),
            self.option("kegg_table_2").prop['path'],
            pathway_file,
            source=self.option('source')
        )
        pngs = os.listdir(pathway_file)
        tar_file = pathway_file + ".tar.gz"
        with tarfile.open(tar_file, mode='w:gz') as f:
            for png in pngs:
                f.add(pathway_file + "/" + png, arcname = png)
        shutil.rmtree(pathway_file)

        # api_geneset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        record_id = self.option("main_table_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            self.error("main_id参数必须为字符串或者ObjectId类型", code = "15000401")
        conn = api_geneset.db["sg_geneset_kegg_class"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            # self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '')
        else:
            # self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '')
        graph_dir = os.path.join(self.workflow_output, 'pathways')
        result_dir = self.workflow_output
        conn.update({"_id": record_id}, {"$set": {'result_dir': result_dir}}, upsert=True)
        self.end()

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if os.path.exists(origin):
            os.link(origin, new)

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        kegg_level = os.path.join(self.work_dir, 'gene_kegg_level_table.xls')
        geneset_kegg_class = os.path.join(self.kegg_class.output_dir, 'kegg_stat.xls')
        chart.prok_geneset_kegg_bar(geneset_kegg_class, kegg_level)
        chart.to_pdf()

        # move pdf
        pdf = glob.glob(os.path.join(self.work_dir, '*.annot_kegg_bar.kegg_bar.pdf'))
        for each in pdf:
            prefix = os.path.basename(each).split('.annot_kegg_bar.kegg_bar.pdf')[0]
            self.move_pdf(each, os.path.join(self.kegg_class.output_dir, prefix + '_kegg_class.pdf'))

    def end(self):
        self.chart()
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "基因集KEGG功能分类结果目录"],
        #     ["kegg_analysis_of_anotate", " ", "KEGG分析结果表"],
        #     ["kegg_statistic", " ", "KEGG统计分析结果表"],
        # ])
        # shutil.rmtree(self.kegg_class.output_dir + "/ko")
        # os.remove(self.kegg_class.output_dir + "/kegg_stat.xls")
        result_dir = self.add_upload_dir(self.kegg_class.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基因集KEGG功能分类结果目录"],
            ["ko", " ", "Pathway代谢通路文件"],
            ["kegg_stat.xls", " ", "Pathway分类统计表"],
            ["pathways.tar.gz", " ", "Pathway代谢通路压缩文件"],
        ])
        result_dir.add_regexp_rules([
            [r".*kegg_class\.pdf", "", "Pathway分类统计柱状图"],
        ])
        super(GenesetKeggWorkflow, self).end()

    def run_kegg_class(self):
        opts = {
            "geneset_kegg": self.option("geneset_kegg").prop['path'],
            "kegg_table": self.option("kegg_table"),
            "geneset_id": self.option("geneset_id"),
            "background_links": self.option("add_info"),
            "type": self.option("type"),
            "task_id": self.option("task_id"),
            "source": self.option('source'),
            "kegg_version": self.option('kegg_version')
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()
