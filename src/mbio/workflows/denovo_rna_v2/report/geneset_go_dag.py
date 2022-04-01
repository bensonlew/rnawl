# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from bson.objectid import ObjectId
import types,os
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
import re
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder


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
        self.go_dag = self.add_tool("denovo_rna_v2.geneset.go_dag")
        self.dump_tool = self.api.api("denovo_rna_v2.api_base")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/08 GO_Dag')
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
        super(GenesetGoDagWorkflow, self).send_log(data)

    def run(self):
        self.go_dag.on("end", self.set_db)
        self.get_run_log()
        self.run_go_dag()
        super(GenesetGoDagWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_geneset_go_dag", main_id=self.option('go_enrich_id'),
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
        self.dump_tool.update_db_record('sg_geneset_go_dag', self.option('go_enrich_id'), output_dir=workflow_output,
                                        status="end", main_id=go_enrich_id)
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.go_dag.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.go_dag.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.go_dag.output_dir, os.path.basename(self.run_log)))
        os.system('cp {} {}'.format(os.path.join(self.go_dag.output_dir,'go_dag.pdf'),self.output_dir))
        os.system('cp {} {}'.format(os.path.join(self.go_dag.output_dir,'go_dag.png'),self.output_dir))
        os.system('cp {} {}'.format(os.path.join(self.go_dag.output_dir,'go_dag.svg'),self.output_dir))
        os.system('cp {} {}'.format(os.path.join(self.go_dag.output_dir,'run_parameter.txt'),self.output_dir))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/08 GO_Dag", "", "GO富集有向无环图", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "GO富集有向无环图文件",0,],
            ["./go_dag.pdf", "", "GO富集有向无环图pdf",0],
            ["./go_dag.png", "", "GO富集有向无环图png",0],
            ["./go_dag.svg", "", "GO富集有向无环图svg",0],
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
