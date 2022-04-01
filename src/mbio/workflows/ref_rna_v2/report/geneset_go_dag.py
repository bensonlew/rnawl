# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import types
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import re
import os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile,linkdir
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
import glob

class GenesetGoDagWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetGoDagWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string", "default": None},
            {"name": "task_id", "type": "string", "default": None},
            {"name": "submit_location", "type": "string", "default": None},
            {"name": "task_type", "type": "string", "default": None},
            {"name": "go_enrich_id", "type": "string", "default": None},
            {"name": "go_enrich_detail", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "go_list", "type": "string", "default": None},
            {"name": "top_num", "type": "string", "default": "20"},
            {"name": "significant_diff", "type": "string", "default": None},
            {"name": "significant_value", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.go_dag = self.add_tool("ref_rna_v2.geneset.go_dag")
        self.dump_tool = self.api.api("ref_rna_v2.go_dag")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 GeneSet/08 GO_Dag')
        self.inter_dirs = []
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="ref_rna_v2",
                                                  main_id=self.option('go_enrich_id'))
            interactiondelete.delete_interactions_records()

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
        super(GenesetGoDagWorkflow, self).send_log(data)

    def run(self):
        self.go_dag.on("end", self.set_db)
        self.get_run_log()
        self.run_go_dag()
        super(GenesetGoDagWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_geneset_go_dag", main_id=self.option('go_enrich_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        workflow_output = self.get_workflow_output_dir()
        go_enrich_id = self.option('go_enrich_id')
        if go_enrich_id is not None:
            if isinstance(go_enrich_id, types.StringTypes):
                go_enrich_id = ObjectId(go_enrich_id)
            elif isinstance(go_enrich_id, ObjectId):
                go_enrich_id = go_enrich_id
        self.dump_tool.run_webroot(main_id=go_enrich_id,
                                   relation_result=self.go_dag.work_dir + "/relation.tsv",
                                   visual=self.go_dag.work_dir + "/relation_visual.json")
        self.dump_tool.update_db_record('sg_geneset_go_dag', self.option('go_enrich_id'), output_dir=workflow_output,
                                        status="end", main_id=go_enrich_id, version="v3.1")
        self.end()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        dag_table = self.go_dag.work_dir + "/relation.tsv"
        chart.chart_geneset_enrich_dag(dag_table)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")[0]
        linkfile(pdf_file, self.go_dag.output_dir + "/go_dag2.pdf")

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.go_dag.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.go_dag.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.go_dag.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.go_dag.output_dir)
        self.inter_dirs = [
            ["03 GeneSet", "", "基因集分析结果目录",0],
            ["03 GeneSet/08 GO_Dag", "", "GO富集有向无环图", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "GO富集有向无环图文件",0,"211542"],
            ["./*.pdf", "", "GO富集有向无环图pdf",0,"211543"],
            ["./go_dag.png", "", "GO富集有向无环图png",0,"211544"],
            ["./go_dag.svg", "", "GO富集有向无环图svg",0,"211545"],
            ["./go_dag.tsv", "", "GO层级关系",0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(GenesetGoDagWorkflow, self).end()

    def run_go_dag(self):
        options = {
            "go_enrich_detail": self.option("go_enrich_detail").prop['path'],
            "go_list": self.option("go_list"),
            "top_num": int(self.option("top_num")),
            "significant_diff": self.option("significant_diff"),
            "significant_value": self.option("significant_value")
        }
        self.go_dag.set_options(options)
        self.go_dag.run()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output
