# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd
import json
import glob
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class ProteinsetCircWorkflow(Workflow):
    """
    蛋白集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetCircWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "enrich_type", "type": "string", "default": "kegg"},
            {"name": "go_type", "type": "string", "deafult":"ALL"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "diff_id", "type": "string"},
            {"name": "enrich_id", "type": "string"},
            {"name": "compare_group", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "diff_fc", "type": "infile", "format": "labelfree.common"},
            {"name": "enrich_table", "type": "infile", "format": "labelfree.common"},
            {"name": "proteinset_id", "type": "string"},
            {"name": "type", "type": "string"} # 指示用origin注释还是latest注释
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.circ = self.add_tool("labelfree.proteinset.enrich2circ")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/04_Enrich/03_Chord')
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
        super(ProteinsetCircWorkflow, self).send_log(data)

    def run(self):
        self.circ.on("end", self.set_db)
        self.run_circ()
        super(ProteinsetCircWorkflow, self).run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        geneset_circ_choose = os.path.join(self.circ.output_dir, "{}_enrich_choose.table".format(self.option('enrich_type').lower()))
        geneset_circ_zscore = os.path.join(self.circ.output_dir,  "enrich_zscore")
        chart.chart_proteinset_circ_web(geneset_circ_choose, geneset_circ_zscore, format(self.option('enrich_type').lower()))
        chart.to_pdf()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_proteinset = self.api.api('dia.proteinset')
        self.logger.info("开始进行富集弦图的导表")
        if self.option("enrich_type") == "GO":
            enrich_circ_file = self.circ.output_dir + '/go_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/go_enrich_detail.table'
        elif self.option("enrich_type") == "KEGG":
            enrich_circ_file = self.circ.output_dir + '/kegg_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/kegg_enrich_detail.table'
        enrich_zscore_file = self.circ.output_dir + '/enrich_zscore'
        api_proteinset.add_circ_graph(self.option("main_table_id"), enrich_circ_file, self.option('enrich_type'))
        api_proteinset.add_circ_detail(self.option("main_table_id"), enrich_detail_file, self.option('enrich_type'))
        api_proteinset.update_circ_main(self.option("main_table_id"), enrich_zscore_file, self.option('enrich_type'))
        self.end()

    def end(self):
        self.chart()
        pdf = glob.glob(os.path.join(self.work_dir, '*pdf'))
        for i in pdf:
            os.link(i, os.path.join(self.circ.output_dir, os.path.basename(i)))

        result_dir = self.add_upload_dir(self.circ.output_dir)
        self.inter_dirs = [
            ["5_Proteinset", "", "蛋白集相关分析结果目录", 0],
            ["5_Proteinset/04_Enrich", "", "蛋白富集相关分析结果目录", 0],
            ["5_Proteinset/04_Enrich/03_Chord", "", "蛋白集弦图结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "蛋白集Circ结果目录"],
            ["./*pdf", "", '富集弦图'],
        ])
        super(ProteinsetCircWorkflow, self).end()

    def run_circ(self):
        opts = {
            "enrich_type" : self.option("enrich_type"),
            "enrich_table" : self.option("enrich_table"),
            "fc_table" : self.option("diff_fc")
        }
        self.circ.set_options(opts)
        self.circ.run()
