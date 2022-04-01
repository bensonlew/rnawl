# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
from collections import defaultdict
import os
import re
import time
from itertools import chain
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.ref_rna_v2.copy_file import CopyFile
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple


class DiffGenesetDoWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffGenesetDoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_do", "type": "string"},
            {"name": "anno_type", "type": "string", "default": "do"},
            {"name": "geneset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "geneset_names", "type": "string"},
            {"name": "type", "type": "string"},
            {'name': 'source', 'type': 'string'},
            {"name": "task_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "do":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Diff_Express/03 DiffExpress_Geneset_Annotion/04 DO/')

        self.do_class = self.add_tool("medical_transcriptome.geneset.do_class")

        self.inter_dirs = []

    def run_do_class(self):
        opts = {
            'do_ids': self.option("geneset_do"),
            'geneset_names': self.option("geneset_names"),
        }
        self.do_class.set_options(opts)
        self.do_class.run()


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
        super(DiffGenesetDoWorkflow, self).send_log(data)

    def run(self):
        self.do_class.on('end', self.set_db)
        self.get_run_log()
        self.run_do_class()
        super(DiffGenesetDoWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_diff_geneset_do_class", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('medical_transcriptome.geneset_annot')
        if self.option("source") == "diff_exp":
            api_geneset = self.api.api('medical_transcriptome.diff_geneset_annot')
        add_name = GeneInfoSupple(id2namedict = None, task_id = self.option("task_id"), level ="G")
        add_columns = [x + " seqs" for x in self.option("geneset_names").split(";")]
        self.logger.info("注释分类分析结果入库")
        if self.option("anno_type") == "do":
            output_file = self.do_class.output_dir + "/do_level2_class.xls"
            try:
                add_name.add_gene_name(output_file, split=";", "do", add_column, self.output_dir + "/do_class_stat.xls")
            except:
                CopyFile().linkfile(output_file, self.output_dir + "/do_class_stat.xls")
            api_geneset.add_geneset_do_detail(output_file, self.option("main_table_id"), geneset_names = self.option("geneset_names") )
        self.end()

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("anno_type") == "do":
            self.inter_dirs = [
                ["01 Diff_Express", "", "差异基因数据挖掘结果目录",0],
                ["01 Diff_Express/03 DiffExpress_Geneset_Annotion", "", "差异基因集功能注释分析", 0],
                ["01 Diff_Express/03 DiffExpress_Geneset_Annotion/04 DO", "", "差异基因集DO功能注释分析", 0]
            ]
            result_dir.add_relpath_rules([
                [".", "", "差异基因集DO功能注释文件", 0],
                ["./do_class_stat.xls", "xls", "DO分类统计表", 0], ],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            )

        super(DiffGenesetDoWorkflow, self).end()


