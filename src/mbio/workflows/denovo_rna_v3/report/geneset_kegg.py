# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import glob
import shutil
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
from mbio.packages.denovo_rna_v2.chart_geneset import ChartGeneset
import glob


class GenesetKeggWorkflow(Workflow):
    """
    基因集功能分类分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetKeggWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_kegg", "type": "string"},
            {"name": "kegg_table", "type": "string"},
            {"name": "kegg_table_2", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "add_info", "type": "string"},  # 底图颜色信息
            {"name": "type", "type": "string"},  # 指示用origin注释还是latest注释
            {"name": "task_id", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
            {"name": "source", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.kegg_class = self.add_tool("denovo_rna_v3.geneset.kegg_class")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/04 KEGG_Annotation')
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
        super(GenesetKeggWorkflow, self).send_log(data)

    def run(self):
        # super(GenesetKeggWorkflow, self).run()
        self.kegg_class.on("end", self.set_db)
        self.get_run_log()
        self.run_kegg_class()
        super(GenesetKeggWorkflow, self).run()
        # self.set_db()
        # self.end()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_geneset_kegg_class", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        # src/mbio/api/database/denovo_rna_v2/geneset.py
        api_geneset = self.api.api('denovo_rna_v3.geneset')
        self.logger.info("开始进行kegg_class的导表")
        output_file = self.kegg_class.output_dir + '/kegg_stat.xls'
        pathway_file = self.kegg_class.output_dir + '/pathways'
        # api_geneset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
        api_geneset.add_kegg_regulate_new(self.option("main_table_id"), self.option("geneset_id"), output_file,
                                          self.option("kegg_table_2"), self.work_dir, self.option("geneset_type"),
                                          source=self.option('source'))
        api_geneset.add_kegg_regulate_pic(self.option("main_table_id"), self.option("kegg_table_2"), pathway_file, source=self.option("source"))
        # api_geneset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        record_id = self.option("main_table_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            self.set_error("main_id参数必须为字符串或者ObjectId类型!", code="13701701")
        conn = api_geneset.db["sg_geneset_kegg_class"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        graph_dir = os.path.join(self.workflow_output, 'pathways')
        conn.update({"_id": record_id}, {"$set": {'graph_dir': graph_dir}}, upsert=True)
        self.end()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"

        kegg_class_table = self.work_dir + "/kegg_statistic"
        chart.chart_geneset_class_kegg(kegg_class_table)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            print "copy", p, self.kegg_class.output_dir + "/" + os.path.basename(p)
            shutil.copyfile(p, self.kegg_class.output_dir + "/" + os.path.basename(p))

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        self.chart()
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "基因集KEGG功能分类结果目录"],
        #     ["kegg_analysis_of_anotate", " ", "KEGG分析结果表"],
        #     ["kegg_statistic", " ", "KEGG统计分析结果表"],
        # ])
        shutil.rmtree(self.kegg_class.output_dir + "/ko")
        os.remove(self.kegg_class.output_dir + "/kegg_stat.xls")
        shutil.copyfile(self.work_dir + "/kegg_analysis_of_anotate",
                        self.kegg_class.output_dir + "/kegg_analysis_detail.xls")
        statis = pd.read_table(self.work_dir + "/kegg_statistic", sep="\t", header=0)
        statis_new = statis.drop(['kegg_id', 'geneset_type', 'geneset_id'], axis=1)
        statis_new.to_csv(self.work_dir + "/kegg_statistics", sep='\t', index=False)
        shutil.copyfile(self.work_dir + "/kegg_statistics", self.kegg_class.output_dir + "/kegg_statistics.xls")

        ## set output
        mark_file = glob.glob(self.kegg_class.output_dir + "/pathways/*.html.mark")
        for file in mark_file:
            os.remove(file)
        os.rename(self.kegg_class.output_dir + '/pathways', self.kegg_class.output_dir + '/kegg_pathways')
        if os.path.exists(os.path.join(self.kegg_class.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.kegg_class.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.kegg_class.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.kegg_class.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/04 KEGG_Annotation", "", "KEGG功能注释", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "KEGG功能注释文件",0,],
            ["kegg_pathways", " ", "KEGG通路图",0],
            ["kegg_analysis_detail.xls", " ", "基因/转录本对应KEGG注释详情表",0],
            ["kegg_statistics.xls", " ", "KEGG分类统计表 ",0,],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ['kegg_annot.*.pdf', 'pdf', 'Pathway分类统计柱状图', 0, ]
        ])
        result_dir.add_regexp_rules([
            [r'kegg_pathways/map.*\.html', '', 'KEGG通路html文件',0,"201476"],
            [r'kegg_pathways/map.*\.png', '', 'KEGG通路png图片',0,"201477"],
        ])
        super(GenesetKeggWorkflow, self).end()

    def run_kegg_class(self):
        opts = {
            "geneset_kegg": self.option("geneset_kegg"),
            "kegg_table": self.option("kegg_table"),
            "kegg_table2": self.option("kegg_table_2"),
            "geneset_id": self.option("geneset_id"),
            "background_links": self.option("add_info"),
            "type": self.option("type"),
            "task_id": self.option("task_id"),
            "kegg_version": self.option('kegg_version'),
            "source": self.option('source')
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()


