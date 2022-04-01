# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import shutil
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json


class GenesetCircWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetCircWorkflow, self).__init__(wsheet_object)
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
            {"name": "diff_fc", "type": "infile", "format": "lnc_rna.common"},
            {"name": "enrich_table", "type": "infile", "format": "lnc_rna.common"},
            {"name": "gene_detail", "type": "string", "default": ""},
            {"name": "anno_list", "type": "string", "default": ""},
            {"name": "geneset_id", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "type", "type": "string"} # 指示用origin注释还是latest注释
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.circ = self.add_tool("lnc_rna.geneset.enrich2circ")
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/04 GeneSet/07 Enrich_Circ')
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
        super(GenesetCircWorkflow, self).send_log(data)

    def run(self):
        self.circ.on("end", self.set_db)
        self.get_run_log()
        self.run_circ()
        super(GenesetCircWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("lnc_rna", table="sg_geneset_circ", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('lnc_rna.lnc_rna_geneset')
        self.logger.info("开始进行富集弦图的导表")
        enrich_circ_file = None
        enrich_detail_file = None
        if self.option("enrich_type") == "GO":
            enrich_circ_file = self.circ.output_dir + '/go_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/go_enrich_detail.table'
        elif self.option("enrich_type") == "KEGG":
            enrich_circ_file = self.circ.output_dir + '/kegg_enrich_choose.table'
            enrich_detail_file = self.circ.output_dir + '/kegg_enrich_detail.table'
        enrich_zscore_file = self.circ.output_dir + '/enrich_zscore'
        if enrich_circ_file:
            api_geneset.add_circ_graph(self.option("main_table_id"), enrich_circ_file, self.option('enrich_type'))
        if enrich_detail_file:
            # geneset_type=None, gene_info=None
            api_geneset.add_circ_detail(self.option("main_table_id"), enrich_detail_file, self.option('enrich_type'),
                                        geneset_type=self.option('geneset_type'),
                                        gene_info=self.option('gene_detail'))
        api_geneset.update_circ_main(self.option("main_table_id"), enrich_zscore_file, self.option('enrich_type'))
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.circ.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.circ.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.circ.output_dir, os.path.basename(self.run_log)))

        if os.path.exists(self.circ.output_dir + '/upload'):
            shutil.rmtree(self.circ.output_dir + '/upload')
        os.mkdir(self.circ.output_dir + '/upload')
        shutil.copyfile(self.circ.output_dir + "/{}_enrich_detail.table".format(self.option("enrich_type").lower()),
                        self.circ.output_dir + "/upload/{}_enrich_detail.table".format(self.option("enrich_type").lower()))
        shutil.copyfile(self.circ.output_dir + "/run_parameter.txt", self.circ.output_dir + "/upload/run_parameter.txt")

        result_dir = self.add_upload_dir(self.circ.output_dir + '/upload')
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录", 0],
            ["04 GeneSet/07 Enrich_Circ", "", "富集弦图", 0],
        ]
        result_dir.add_relpath_rules([
            [".", "", "富集弦图文件", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'.*_enrich_detail\.table', 'table', '富集弦图数据表', 0],
        ])
        super(GenesetCircWorkflow, self).end()

    def run_circ(self):
        opts = {
            "enrich_type" : self.option("enrich_type"),
            "enrich_table" : self.option("enrich_table"),
            "fc_table" : self.option("diff_fc"),
            "p_thre" : self.option("p_thre"),
            "padj_thre" : self.option("padj_thre"),
            "anno_num_thre" : self.option("anno_num_thre"),
            "anno_list" : self.option("anno_list"),
            "gene_detail" : self.option("gene_detail")
        }
        self.circ.set_options(opts)
        self.circ.run()
