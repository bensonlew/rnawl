# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
from biocluster.config import Config
import glob
import os
import re
from bson.objectid import ObjectId
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import shutil


class GenesetEnrichWorkflow(Workflow):
    """
    基因集富集分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(GenesetEnrichWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "geneset_kegg", "type": "string"},
            {"name": "kegg_table", "type": "string"},
            {"name": "go_list", "type": "string"},
            {'name': 'go_version', 'type': 'string', 'default': '2019'},
            {"name": "geneset_list", "type": "string"},
            {"name": "all_list", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "type", "type": "string"},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
            {"name": "add_info", "type": "string", "default": None},  # 输入两列的列表文件，有head，第一列为pathway，第二列为底图链接
            {"name": "geneset_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.enrich_tool = self.add_tool("lnc_rna.geneset.go_enrich") \
            if self.option("anno_type") == "go" else self.add_tool("lnc_rna.geneset.kegg_rich")
        self.kegg_class = self.add_tool("lnc_rna.geneset.kegg_class")
        self.output_dir1 = self.enrich_tool.output_dir
        self.work_dir2 = self.kegg_class.work_dir
        if self.option("anno_type") == "go":
            self._sheet.output = self._sheet.output.replace('interaction_results',
                                                            'interaction_results/04 GeneSet/05 GO_Enrich')
        elif self.option("anno_type") == "kegg":
            self._sheet.output = self._sheet.output.replace('interaction_results',
                                                            'interaction_results/04 GeneSet/06 KEGG_Enrich')
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
        super(GenesetEnrichWorkflow, self).send_log(data)

    def run(self):
        if self.option("anno_type") == "kegg":
            options = {
                "kegg_table": self.option("kegg_table"),
                # "all_list": background_path,
                "diff_list": self.option("geneset_list"),
                "correct": self.option("method"),
                "add_info": self.option("add_info")
            }
        else:
            options = {
                "diff_list": self.option("geneset_list"),
                # "all_list": background_path,
                'go_version': self.option('go_version'),
                "go_list": self.option("go_list"),
                # "pval": self.option("pval"),
                "method": self.option("method"),
            }
        self.logger.info(options)
        self.enrich_tool.set_options(options)
        if self.option("anno_type") == "kegg":
            self.enrich_tool.on('end', self.run_kegg_class)
            self.kegg_class.on('end', self.set_db)
        else:
            self.enrich_tool.on('end', self.set_db)
        self.get_run_log()
        self.enrich_tool.run()
        super(GenesetEnrichWorkflow, self).run()

    def get_run_log(self):
        if self.option("anno_type") == "go":
            table = "sg_geneset_go_enrich"
        else:
            table = "sg_geneset_kegg_enrich"
        get_run_log = GetRunLog("lnc_rna", table=table, main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('lnc_rna.lnc_rna_geneset')
        output_file = glob.glob("{}/*.xls".format(self.enrich_tool.output_dir))[0]
        workflow_output = self.get_workflow_output_dir() + '/' + self.option("anno_type") + '_enrich_stat.xls'
        # png_file = glob.glob("{}/*.png".format(self.output_dir))
        # go_png = self.output_dir + "/go_lineage.png"
        # go_pdf = self.output_dir + "/go_lineage.pdf"
        go_adjust_png = self.output_dir1 + "/adjust_lineage.png"
        go_adjust_pdf = self.output_dir1 + "/adjust_lineage.pdf"
        if self.option("anno_type") == "kegg":

            api_geneset.add_kegg_enrich_detail(self.option("main_table_id"), output_file,
                                               self.option("geneset_list"), self.option("all_list"))
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option("main_table_id"),
                                         result_dir=workflow_output)
            api_geneset.add_kegg_enrich_pic(self.option("main_table_id"), output_file,
                                            self.kegg_class.output_dir + "/pathways")
            graph_dir = os.path.join(self.get_workflow_output_dir(), 'pathways')
            api_geneset.update_db_record('sg_geneset_kegg_enrich', self.option("main_table_id"),
                                         graph_dir=graph_dir, status="end")
        else:
            api_geneset.add_go_enrich_detail(self.option("main_table_id"), output_file)
            # api_geneset.update_directed_graph(self.option("main_table_id"), go_adjust_png, go_adjust_pdf)
            api_geneset.update_db_record('sg_geneset_go_enrich', self.option("main_table_id"),
                                         result_dir=workflow_output, status="end")
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))

        if self.option("anno_type") == "go":
            os.rename(self.output_dir1 + "/go_enrich_geneset_list_gene.xls", self.output_dir + "/go_enrich_stat.xls")
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录", 0],
                ["04 GeneSet/05 GO_Enrich", "", "GO功能富集", 0],
            ]
            result_dir.add_relpath_rules([
                [".", "", "GO富集分析文件", 0, "211073"],
                ["./go_enrich_stat.xls", "xls", "GO富集分析统计表"],
                ['run_parameter.txt', 'txt', '运行参数日志', 0]
            ])
        elif self.option("anno_type") == "kegg":
            kegg_file = glob.glob(self.output_dir1 + "/*kegg*.xls")
            shutil.copyfile(kegg_file[0], self.output_dir + "/kegg_enrich_stat.xls")
            shutil.copytree(self.work_dir2 + "/png", self.output_dir + "/kegg_enrich_pathways")
            result_dir = self.add_upload_dir(self.output_dir)
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录", 0],
                ["04 GeneSet/06 KEGG_Enrich", "", "KEGG富集分析", 0],
            ]
            result_dir.add_relpath_rules([
                [".", "", "KEGG富集分析文件", 0, "211074"],
                ["./kegg_enrich_pathways", " ", "KEGG富集通路图", 0, "211075"],
                ["./kegg_enrich_stat.xls", "xls", "KEGG富集分析统计表"],
                ['run_parameter.txt', 'txt', '运行参数日志', 0]
            ])
            result_dir.add_regexp_rules([
                [r'.*/map.*\.png', 'png', 'KEGG通路图片png', 0],
                [r'.*/map.*\.html', 'html', 'KEGG通路图片html', 0],
            ])

        super(GenesetEnrichWorkflow, self).end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:', '/mnt/ilustre/data/')
        return workflow_output

    def run_kegg_class(self):
        opts = {
            "geneset_kegg": self.option("geneset_kegg"),
            "kegg_table": self.option("kegg_table"),
            "geneset_id": self.option("geneset_id"),
            "background_links": self.option("add_info"),
            "type": self.option("type"),
            "task_id": self.option("task_id"),
            "kegg_version": self.option('kegg_version')
        }
        self.kegg_class.set_options(opts)
        self.kegg_class.run()
