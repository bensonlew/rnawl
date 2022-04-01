#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/30 10:25
@file    : geneset_gsea.py
"""

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import re
from bson.objectid import ObjectId
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class GenesetGseaWorkflow(Workflow):
    """
    基因集富集分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GenesetGseaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "task_id", "type": "string"},

            {"name": "preranked", "type": "bool", 'default': False},
            {"name": "geneset_source", "type": "string"},
            {"name": "geneset_id", "type": "string"},

            {"name": "g1", "type": "string"},
            {"name": "g2", "type": "string"},
            {"name": "group", "type": "string"},
            {"name": "group_dict", "type": "string"},
            {"name": "group_id", "type": "string"},

            {"name": "max_num", "type": "string", "default": "500"},
            {"name": "min_num", "type": "string", "default": "15"},

            {"name": "level", "type": "string", "default": "G"},
            {"name": "go_genesets", "type": "string", "default": "all"},
            {"name": "go_list", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "go_type", "type": "string", "default": ""},

            {"name": "sort_method", "type": "string"},
            {"name": "permutation_type", "type": "string"},
            {"name": "plot_top_x", "type": "string", "default": "20"},



            {"name": "gmx", "type": "infile", "format": "ref_rna_v2.geneset_gmt"},
            {"name": "rnk", "type": "infile", "format": "ref_rna_v2.geneset_rnk"},
            {"name": "matrix", "type": "infile", "format": "ref_rna_v2.ref_common"},
            {"name": "genes_detail", "type": "infile", "format": "ref_rna_v2.ref_common"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.step.add_steps('gsea_analysis')
        if self.option("go_list").is_set:
            self.go_gene_gmx = self.add_tool('ref_rna_v2.geneset.go_gene_gmx')
        if self.option('preranked'):
            self.gsea_tool = self.add_tool('ref_rna_v2.geneset.gsea_preranked')
        else:
            self.gsea_tool = self.add_tool('ref_rna_v2.geneset.gsea')
        self.gsea_tool.on('end', self.set_db)
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/07 Advanced_Analysis/06 GSEA_Analysis')
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
        super(GenesetGseaWorkflow, self).send_log(data)

    def get_go_gmx(self):
        opts = {
            "go_list": self.option("go_list").prop['path'],
            "min_num": self.option("min_num"),
            "max_num": self.option("max_num"),
            "go_type": self.option("go_type"),
            "go_sets": self.option("go_genesets")
        }
        self.go_gene_gmx.set_options(opts)
        self.go_gene_gmx.run()



    def gsea_runner(self):
        # 输出: biomart.xls, biomart.json
        if self.option("go_list").is_set:
            gmx = self.go_gene_gmx.output_dir + '/go.gmt'
        else:
            gmx =self.option('gmx')
        self.gsea_tool.set_options({
            'g1': self.option('g1'),
            'g2': self.option('g2'),
            'group': self.option('group'),
            'set_max': int(self.option('max_num')),
            'set_min': int(self.option('min_num')),
            'metric': self.option('sort_method'),
            'plot_top_x': int(self.option('plot_top_x')),
            'permute': self.option('permutation_type'),
            'geneset_source': self.option('geneset_source'),
            'gmx': gmx,
            'matrix': self.option('matrix'),
            'gene_detail': self.option('genes_detail')
        })
        self.gsea_tool.on('start', self.set_step, {'start': self.step.gsea_analysis})
        self.gsea_tool.on('end', self.set_step, {'end': self.step.gsea_analysis})
        self.gsea_tool.run()

    def gsea_preranked_runner(self):
        if self.option("go_list").is_set:
            gmx = self.go_gene_gmx.output_dir + '/go.gmt'
        else:
            gmx =self.option('gmx')
        self.gsea_tool.set_options({
            'set_max': self.option('max_num'),
            'set_min': self.option('min_num'),
            'gmx': gmx,
            'rnk': self.option('rnk'),
            'geneset_source': self.option('geneset_source'),
            'matrix': self.option('matrix'),
            'plot_top_x': int(self.option('plot_top_x')),
            'gene_detail': self.option('genes_detail')
        })
        self.gsea_tool.on('start', self.set_step, {'start': self.step.gsea_analysis})
        self.gsea_tool.on('end', self.set_step, {'end': self.step.gsea_analysis})
        self.gsea_tool.run()

    def run(self):
        self.get_run_log()
        if self.option("preranked"):
            if self.option("go_list").is_set:
                self.go_gene_gmx.on('end', self.gsea_preranked_runner)
                self.get_go_gmx()
            else:
                self.gsea_preranked_runner()
        else:
            if self.option("go_list").is_set:
                self.go_gene_gmx.on('end', self.gsea_runner)
                self.get_go_gmx()
            else:
                self.gsea_runner()
        super(GenesetGseaWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="geneset_gsea", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('whole_transcriptome.gsea')
        api_geneset.run_webroot(self.option("main_table_id"), self.gsea_tool.output_dir)

        conn = api_geneset.db["geneset_gsea"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')

        graph_dir = self.workflow_output
        conn.update({"_id": ObjectId(self.option("main_table_id"))}, {"$set": {'result_dir': graph_dir}}, upsert=True)

        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.gsea_tool.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.gsea_tool.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.gsea_tool.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.gsea_tool.output_dir)
        self.inter_dirs = [
            ["07 Advanced_Analysis", "", "高级分析结果目录",0],
            ["07 Advanced_Analysis/06 GSEA_Analysis", "", "GSEA分析结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "GSEA分析文件",0,"211546"],
            ["./*.pdf", "", "该基因集富集结果图pdf",0,"211547"],
            ["./*.png", "", "该基因集富集结果图png",0,"211548"],
            ["./all_sets.detail", "txt", "所有基因集富集结果详情表",0,"211549"],
            ["./gsea_report.xls", "xls", "所有基因集富集结果统计表",0,"211550"],
            ["./all_exp.detail", "txt", "所有基因集leading基因表达量信息",0,"211551"],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(GenesetGseaWorkflow, self).end()
