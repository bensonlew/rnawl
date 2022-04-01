# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd
import glob
import json
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder

class ProteinsetKeggWorkflow(Workflow):
    """
    蛋白集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetKeggWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "proteinset_kegg", "type": "infile", "format": "labelfree.common"},
            {"name": "kegg_table", "type": "infile", "format": "labelfree.kegg_table"},
            {"name": "kegg_table_2", "type": "infile", "format": "denovo_rna_v2.vcf"},
            {"name": "proteinset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "proteinset_id", "type": "string"},
            {"name": "add_info", "type": "string"},  # 底图颜色信息
            {"name": "type", "type": "string"}, # 指示用origin注释还是latest注释
            {"name": "task_id", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': "2017"},
            {'name': 'task_version', 'type': 'string', 'default': "2"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.kegg_class = self.add_tool("dia_v3.proteinset.kegg_class")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/5_Proteinset/03_Anno/02_KEGG')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(ProteinsetKeggWorkflow, self).send_log(data)

    def run(self):
        # super(ProteinsetKeggWorkflow, self).run()
        self.kegg_class.on("end", self.set_db)
        self.run_kegg_class()
        super(ProteinsetKeggWorkflow, self).run()
        # self.set_db()
        # self.end()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_proteinset = self.api.api('dia.proteinset')

        kegg_stat = self.kegg_class.output_dir + '/kegg_stat.xls'
        pathway_class = self.option("kegg_table_2").prop['path']
        pathway_file = self.kegg_class.output_dir + '/pathways'

        self.logger.info("开始进行kegg_class的导表")
        # api_proteinset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
        #api_proteinset.add_kegg_regulate_new(self.option("main_table_id"), self.option("proteinset_id"), output_file, self.option("kegg_table_2").prop['path'], self.work_dir)
        api_proteinset.add_kegg_regulate_new2(self.option("main_table_id"), self.option("proteinset_kegg").prop['path'], kegg_stat, self.option("kegg_table_2").prop['path'])
        api_proteinset.add_kegg_regulate_pic(self.option("main_table_id"), self.option("kegg_table_2").prop['path'],
                                             pathway_file)

        # api_proteinset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        record_id = self.option("main_table_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = api_proteinset.db["sg_proteinset_kegg_class"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        graph_dir = os.path.join(self.workflow_output, 'pathways')
        result_dir = self.workflow_output
        conn.update({"_id": record_id}, {"$set": {'result_dir': result_dir}}, upsert=True)
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        kegg_class_stat_updown = os.path.join(self.kegg_class.output_dir, "kegg_stat.xls")
        kegg_class_level = os.path.join(self.work_dir, "protein_kegg_level_table.xls")
        chart.chart_proteinsetcluster_keggclass_web(kegg_class_stat_updown, kegg_class_level)
        chart.to_pdf()

        # move pdf to result dir
        pdf_files = glob.glob(self.work_dir + "/*kegg.barline.pdf")
        for pdf_file in pdf_files:
            os.link(pdf_file, os.path.join(self.kegg_class.output_dir, os.path.basename(pdf_file).split('.barline')[0]+'.pdf'))

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "蛋白集KEGG功能分类结果目录"],
        #     ["kegg_analysis_of_anotate", " ", "KEGG分析结果表"],
        #     ["kegg_statistic", " ", "KEGG统计分析结果表"],
        # ])
        # shutil.rmtree(self.kegg_class.output_dir + "/ko")
        # os.remove(self.kegg_class.output_dir + "/kegg_stat.xls")
        self.chart()
        rm_files = glob.glob(os.path.join(self.kegg_class.output_dir, 'pathways/*.png')) + \
                   glob.glob(os.path.join(self.kegg_class.output_dir, 'pathways/*.pdf'))
        if self.option("task_version") >= "2.1":
            for rm_file in rm_files:
                os.remove(rm_file)
        result_dir = self.add_upload_dir(self.kegg_class.output_dir)
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集分析", 0],
            ["5_Proteinset/03_Anno/", "", "功能注释", 0],
            ["5_Proteinset/03_Anno/02_KEGG", "", "KEGG注释", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白集KEGG功能分类结果目录", 0, "220079"],
            ["pathways", " ", "KEGG分析结果通路表", 0, "220080"],
            ["kegg_analysis_of_anotate.xls", " ", "KEGG分类统计表", 0, "220081"],
            ["kegg_statistics.xls", " ", "KEGG分类统计结果表", 0, "220082"],
            ["*kegg.pdf", " ", "Pathway分类统计柱状图", 0, "220083"],
        ])
        super(ProteinsetKeggWorkflow, self).end()

    def run_kegg_class(self):
        opts = {
            "proteinset_kegg": self.option("proteinset_kegg").prop['path'],
            "kegg_table": self.option("kegg_table"),
            "proteinset_id": self.option("proteinset_id"),
            "background_links": self.option("add_info"),
            "type": self.option("type"),
            "task_id": self.option("task_id"),
            "kegg_version": self.option('kegg_version')
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()
