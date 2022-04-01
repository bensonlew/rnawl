# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import glob
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json


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
            {'name': 'kegg_version', 'type': 'string', 'default': "2017"},
            {"name": "source", "type": "string", "default": ""}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.kegg_class = self.add_tool("lnc_rna.geneset.kegg_class")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/04 GeneSet/04 KEGG_Annotation')
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
        self.kegg_class.on("end", self.set_db)
        self.get_run_log()
        self.run_kegg_class()
        super(GenesetKeggWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_geneset_kegg_class", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('lnc_rna.lnc_rna_geneset')
        self.logger.info("开始进行kegg_class的导表")
        output_file = self.kegg_class.output_dir + '/kegg_stat.xls'
        pathway_file = self.kegg_class.output_dir + '/pathways'
        # api_geneset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
        api_geneset.add_kegg_regulate_new(self.option("main_table_id"), self.option("geneset_id"), output_file,
                                          self.option("kegg_table_2"), self.work_dir, self.option("geneset_type"))
        api_geneset.add_kegg_regulate_pic(self.option("main_table_id"), self.option("kegg_table_2"), pathway_file)
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
        # if re.match(r'tsanger:', self.workflow_output_tmp):
        #     self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        # else:
        #     self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        # graph_dir = os.path.join(self.workflow_output, 'pathways')
        # 去除更新 'graph_dir': graph_dir 状态的字段
        conn.update({"_id": record_id}, {"$set": {'status': 'end'}}, upsert=True)
        self.end()

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "基因集KEGG功能分类结果目录"],
        #     ["kegg_analysis_of_anotate", " ", "KEGG分析结果表"],
        #     ["kegg_statistic", " ", "KEGG统计分析结果表"],
        # ])
        shutil.rmtree(self.kegg_class.output_dir + "/ko")
        os.remove(self.kegg_class.output_dir + "/kegg_stat.xls")
        shutil.copyfile(self.work_dir + "/kegg_analysis_of_anotate",
                        self.kegg_class.output_dir + "/pathway_class_stat.xls")
        statis = pd.read_table(self.work_dir + "/kegg_statistic", sep="\t", header=0)
        statis_new = statis.drop(['kegg_id', 'geneset_type', 'geneset_id'], axis=1)
        statis_new.to_csv(self.work_dir + "/kegg_statistics", sep='\t', index=False)
        shutil.copyfile(self.work_dir + "/kegg_statistics", self.kegg_class.output_dir + "kegg_class_stat.xls")
        if os.path.exists(os.path.join(self.kegg_class.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.kegg_class.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.kegg_class.output_dir, os.path.basename(self.run_log)))

        shutil.copyfile(self.kegg_class.output_dir + '/pathway_class_stat.xls', self.output_dir + '/pathway_class_stat.xls')
        shutil.copyfile(self.kegg_class.output_dir + '/run_parameter.txt', self.output_dir + '/run_parameter.txt')
        shutil.copytree(self.kegg_class.work_dir + '/png', self.output_dir + '/pathways')

        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录", 0],
            ["04 GeneSet/04 KEGG_Annotation", "", "KEGG功能注释", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "KEGG功能注释文件", 0],
            ['./kegg_class_stat.xls', 'xls', 'KEGG分类统计表', 0],
            ["./pathway_class_stat.xls", "xls", "Pathway分类统计表",0],
            ["./pathways", "", "KEGG通路图", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'.*/map.*\.png', 'png', 'KEGG通路图片png', 0],
            [r'.*/map.*\.html', 'html', 'KEGG通路图片html', 0],
        ])
        super(GenesetKeggWorkflow, self).end()

    def run_kegg_class(self):
        opts = {
            "geneset_kegg": self.option("geneset_kegg"),
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
