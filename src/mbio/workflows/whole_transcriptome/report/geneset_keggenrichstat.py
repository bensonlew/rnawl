# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd
from bson.objectid import ObjectId
from collections import defaultdict
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class GenesetKeggenrichstatWorkflow(Workflow):
    """
    基因集功能分类分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetKeggenrichstatWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_kegg_enrich_info", "type": "string"},
            {"name": "geneset_kegg_enrich_stat_id", "type": "string"},
            #{"name": "draw_type", "type": "string"},
            {"name":"pathway_ids","type":"string"},
            {"name": "update_info", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "stat_level", "type": "string","default":"pvalue"},
            {"name": "stat_threshold_value", "type": "float","defalut":0.05},
            {"name": "stat_numbers_value", "type": "int","default":20},
            #{"name": "geneset_name", "type": "string"},
            {"name": "geneset_list", "type": "string"},  # 底图颜色信息
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("whole_transcriptome.geneset.kegg_richstat")

    def run(self):
        # super(GenesetKeggWorkflow, self).run()
        self.tool.on("end", self.set_db)
        self.get_run_log()
        self.run_tool()
        super(GenesetKeggenrichstatWorkflow,self).run()
        # self.set_db()
        # self.end()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="geneset_kegg_enrich_stat",
                                main_id=self.option('geneset_kegg_enrich_stat_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('whole_transcriptome.whole_transcriptome_geneset')
        self.logger.info("开始进行kegg_enrich_class的导表")
        output_file = self.tool.work_dir + "/" + "kegg_statistic"

        api_geneset.add_genesetkegg_enrich_stat(self.option("geneset_kegg_enrich_stat_id"),output_file)
        record_id = self.option("geneset_kegg_enrich_stat_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            self.set_error("main_id参数必须为字符串或者ObjectId类型!", code="13701701")
        self.end()

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "基因集KEGG功能分类结果目录"],
        #     ["kegg_analysis_of_anotate", " ", "KEGG分析结果表"],
        #     ["kegg_statistic", " ", "KEGG统计分析结果表"],
        # ])
        statis = pd.read_table(self.tool.work_dir + "/kegg_statistic", sep="\t", header=0)
        shutil.copyfile(self.tool.work_dir + "/kegg_statistic", self.output_dir + "/kegg_enrich_statistic.xls")
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "基因集KEGG富集功能分类结果目录", 0, "211080"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(GenesetKeggenrichstatWorkflow, self).end()


    def run_tool(self):
        opts = {
            "geneset_kegg_enrich_info": self.option("geneset_kegg_enrich_info"),
            "geneset_kegg_enrich_stat_id": self.option("geneset_kegg_enrich_stat_id"),
            #"geneset_name": self.option("geneset_name"),
        }
        self.tool.set_options(opts)
        self.tool.run()

