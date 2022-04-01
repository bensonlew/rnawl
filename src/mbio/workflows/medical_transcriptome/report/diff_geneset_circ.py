# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import glob
import shutil
import json
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset


class DiffGenesetCircWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffGenesetCircWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "enrich_type", "type": "string", "default": "kegg"},
            {"name": "go_type", "type": "string", "deafult":"ALL"},
            {"name": "main_table_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "diff_id", "type": "string"},
            {"name": "p_thre", "type": "float"},
            {"name": "padj_thre", "type": "float"},
            {"name": "anno_num_thre", "type": "int"},
            {"name": "enrich_id", "type": "string"},
            {"name": "compare_group", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "diff_fc", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "enrich_table", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "gene_detail", "type": "string", "default": ""},
            # 2019.01.18 bug anno_list的type为string而不是int
            {"name": "anno_list", "type": "string", "default": ""},
            {"name": "geneset_id", "type": "string"},
            {"name": "type", "type": "string"} # 指示用origin注释还是latest注释
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.circ = self.add_tool("medical_transcriptome.geneset.geneset_enrich2circ")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Diff_Express/04 DiffExpress_Geneset_Enrich/06 Enrich_Circ')

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
        super(DiffGenesetCircWorkflow, self).send_log(data)

    def run(self):
        self.circ.on("end", self.set_db)
        self.get_run_log()
        self.run_circ()
        super(DiffGenesetCircWorkflow, self).run()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        if self.option("enrich_type").upper() == "GO":
            circ_table = self.circ.output_dir + '/go_enrich_choose.table'
        elif self.option("enrich_type").upper() == "REACTOME":
            circ_table = self.circ.output_dir + '/reactome_enrich_choose.table'
        elif self.option("enrich_type").upper() == "DO":
            circ_table = self.circ.output_dir + '/do_enrich_choose.table'
        elif self.option("enrich_type").upper() == "DISGENET":
            circ_table = self.circ.output_dir + '/disgenet_enrich_choose.table'
        else:
            circ_table = self.circ.output_dir + '/kegg_enrich_choose.table'
        circ_zscore_input = self.circ.output_dir + '/enrich_zscore'
        chart.chart_geneset_enrich_circ(circ_table, circ_zscore_input)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")[0]
        os.link(pdf_file, self.circ.output_dir + "/enrich_circ.pdf")

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_diff_geneset_circ", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('medical_transcriptome.diff_geneset')
        self.logger.info("开始进行富集弦图的导表")
        if self.option("enrich_type") == "GO":
            enrich_circ_file = self.circ.output_dir + '/go_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/go_enrich_detail.table'
        elif self.option("enrich_type") == "KEGG":
            enrich_circ_file = self.circ.output_dir + '/kegg_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/kegg_enrich_detail.table'
        elif self.option("enrich_type").upper() == "DISGENET":
            enrich_circ_file = self.circ.output_dir + '/disgenet_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/disgenet_enrich_detail.table'
        elif self.option("enrich_type").upper() == "DO":
            enrich_circ_file = self.circ.output_dir + '/do_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/do_enrich_detail.table'
        elif self.option("enrich_type").upper() == "REACTOME":
            enrich_circ_file = self.circ.output_dir + '/reactome_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/reactome_enrich_detail.table'
        enrich_zscore_file = self.circ.output_dir + '/enrich_zscore'
        api_geneset.add_circ_graph(self.option("main_table_id"), enrich_circ_file, self.option('enrich_type'))
        api_geneset.add_circ_detail(self.option("main_table_id"), enrich_detail_file, self.option('enrich_type'))
        api_geneset.update_circ_main(self.option("main_table_id"), enrich_zscore_file, self.option('enrich_type'))
        self.end()

    def end(self):
        self.chart()
        zscore_file = glob.glob(self.circ.output_dir + "/enrich_zscore")
        go_file = glob.glob(self.circ.output_dir + "/go_enrich_choose.table")
        kegg_file = glob.glob(self.circ.output_dir + "/kegg_enrich_choose.table")
        do_file = glob.glob(self.circ.output_dir + "/do_enrich_choose.table")
        reactome_file = glob.glob(self.circ.output_dir + "/reactome_enrich_choose.table")
        disgenet_file = glob.glob(self.circ.output_dir + "/disgenet_enrich_choose.table")
        rm_file = zscore_file + go_file + kegg_file + do_file + reactome_file + disgenet_file
        for file in rm_file:
            os.remove(file)
        if os.path.exists(os.path.join(self.circ.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.circ.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.circ.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.circ.output_dir)
        self.inter_dirs = [
            ["01 Diff_Express", "", "差异基因数据挖掘结果目录", 0],
            ["01 Diff_Express/04 DiffExpress_Geneset_Enrich", "", "差异基因集功能富集分析", 0],
            ["01 Diff_Express/04 DiffExpress_Geneset_Enrich/06 Enrich_Circ", "", "差异基因富集弦图", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "富集弦图文件",0],
            ["go_enrich_detail.table", "", "富集弦图数据表",0],
            ['./run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'.*\.pdf', '', '富集弦图', 0]
        ])
        super(DiffGenesetCircWorkflow, self).end()

    def run_circ(self):
        opts = {
            "enrich_type" : self.option("enrich_type"),
            "enrich_table" : self.option("enrich_table").prop["path"],
            "fc_table" : self.option("diff_fc").prop["path"],
            "p_thre" : self.option("p_thre"),
            "padj_thre" : self.option("padj_thre"),
            "anno_num_thre" : self.option("anno_num_thre"),
            "anno_list" : self.option("anno_list"),
            "gene_detail" : self.option("gene_detail")
        }
        self.circ.set_options(opts)
        self.circ.run()
