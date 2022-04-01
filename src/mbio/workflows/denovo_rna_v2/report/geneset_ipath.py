# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import json
import shutil
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class GenesetIpathWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetIpathWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "geneset_kegg", "type": "string"},
            {"name": "kegg_table", "type": "infile", "format": "ref_rna_v2.kegg_table"},
            {"name": "main_table_id", "type": "string"},
            {"name": "t2g", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "type", "type": "string" , "default": "origin"} # 指示用origin注释还是latest注释
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ipath = self.add_tool("lnc_rna.geneset.ipath")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/06 Advanced_Analysis/05 iPath_Analysis')
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
        super(GenesetIpathWorkflow, self).send_log(data)

    def run(self):
        self.ipath.on("end", self.set_db)
        self.get_run_log()
        self.run_ipath()
        super(GenesetIpathWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_geneset_ipath", main_id=self.option('main_table_id'), dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('denovo_rna_v2.geneset')
        self.logger.info("开始进行kegg_class的导表")
        output_file = self.ipath.output_dir + '/gene_ipath_input.xls'

        if self.option("geneset_type") == "T":
            api_geneset.add_ipath_detail(self.option("main_table_id"), output_file, self.option('geneset_kegg'), self.option("t2g").prop['path'])
        else:
            api_geneset.add_ipath_detail(self.option("main_table_id"), output_file, self.option('geneset_kegg'))
        # api_geneset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        record_id = self.option("main_table_id")
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
            record_id = record_id
        else:
            self.set_error("main_id参数必须为字符串或者ObjectId类型!", code = "12001401")
        conn = api_geneset.db["sg_geneset_ipath"]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
        graph_dir = self.workflow_output
        conn.update({"_id": record_id}, {"$set": {'result_dir': graph_dir}}, upsert=True)
        self.end()

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        if os.path.exists(os.path.join(self.ipath.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.ipath.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.ipath.output_dir, os.path.basename(self.run_log)))
        os.system('cp {} {}'.format(os.path.join(self.ipath.output_dir,'gene_ipath_input.xls'),self.output_dir))
        os.system('cp {} {}'.format(os.path.join(self.ipath.output_dir,'run_parameter.txt'),self.output_dir))
        result_dir = self.add_upload_dir(self.ipath.output_dir)
        self.inter_dirs = [
            ["06 Advanced_Analysis", "", "高级分析结果目录",0],
            ["06 Advanced_Analysis/05 iPath_Analysis", "", "iPath代谢通路分析结果目录", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "iPath代谢通路分析文件", 0],
            ["./gene_ipath_input.xls", "", "代谢通路数据",0,"211552"],
            # ["./Secondary_metabolites.svg", "", "生物合成次生代谢产物的途径",0],
            # ["./Microbial_metabolism.svg", "", "微生物在不同环境中的代谢", 0],
            # ["./Metabolism.svg", "", "生物系统完整的代谢系统", 0],
            # ["./Antibiotics.svg", "", "各种抗生素的生物合成途径", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(GenesetIpathWorkflow, self).end()

    def run_ipath(self):
        opts = {
            "geneset_kegg": self.option("geneset_kegg"),
            "kegg_table": self.option("kegg_table").prop['path'],
        }
        self.ipath.set_options(opts)
        self.ipath.run()
