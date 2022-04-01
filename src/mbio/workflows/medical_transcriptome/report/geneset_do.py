# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
from collections import defaultdict
import os
import re
import time
import glob
from itertools import chain
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.ref_rna_v2.copy_file import CopyFile
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile

class GenesetDoWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetDoWorkflow, self).__init__(wsheet_object)
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
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/02 Annotation/04 DO')

        self.do_class = self.add_tool("medical_transcriptome.geneset.do_class")
        self.has_results = False
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="medical_transcriptome",
                                                  main_id=self.option('main_table_id'))
            interactiondelete.delete_interactions_records()
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
        super(GenesetDoWorkflow, self).send_log(data)

    def run(self):
        self.do_class.on('end', self.set_db)
        self.get_run_log()
        self.run_do_class()
        super(GenesetDoWorkflow, self).run()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        output_file = self.do_class.output_dir + "/do_level2_class.xls"
        if self.has_results:
            chart.chart_geneset_class_do(do_class_table=output_file)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = glob.glob(self.work_dir + "/*.pdf")
            for p in pdf_file:
                linkfile(p, self.output_dir + "/" + os.path.basename(p))
        else:
            pass

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_geneset_do_class", main_id=self.option('main_table_id'),
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
        add_columns = [x + " seqs" for x in self.option("geneset_names").split(",")]

        self.logger.info("注释分类分析结果入库")
        output_file = self.do_class.output_dir + "/do_level2_class.xls"
        results_num = len(open(output_file, 'r').readlines())
        if results_num > 1:
            self.has_results = True
            if self.option("anno_type") == "do":
                try:
                    add_name.add_gene_name(output_file, split=";", annot_type="do", add_columns=add_columns, outfile_path=self.output_dir + "/do_class_stat.xls")
                except Exception as e:
                    self.logger.info("添加基因name 失败 {}".format(e))
        # CopyFile().linkfile(output_file, self.output_dir + "/do_class_stat.xls")
            api_geneset.add_geneset_do_detail(output_file, self.option("main_table_id"), geneset_names = self.option("geneset_names") )
        else:
            api_geneset.update_db_record('sg_geneset_do_class', self.option("main_table_id"), has_results=False,
                                  status="end")
        self.end()

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("anno_type") == "do":
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录",0],
                ["04 GeneSet/02 Annotation", "", "基因集功能注释", 0],
                ["04 GeneSet/02 Annotation/04 DO", "", "基因集DO功能注释", 0]
            ]
            result_dir.add_relpath_rules([
                [".", "", "基因集DO功能注释文件",0],
                ["*.pdf", "pdf", "DO分类条形图", 0],
                ["./do_class_stat.xls", "xls", "基因集DO分类统计表", 0],
                ['./run_parameter.txt', 'txt', '运行参数日志', 0]
        ])

        super(GenesetDoWorkflow, self).end()


