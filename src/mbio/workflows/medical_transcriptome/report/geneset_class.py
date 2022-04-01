# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

from biocluster.workflow import Workflow
from collections import defaultdict
import os
import re
import time
from itertools import chain
# from mbio.packages.medical_transcriptome.kegg_regulate import KeggRegulate
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset
import glob
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile

class GenesetClassWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetClassWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_go", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            # {"name": "type", "type": "string"},
            {"name": "task_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "go":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/02 Annotation/01 GO')

        self.inter_dirs = []
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete =  InteractionDelete(bind_object=self,project_type="medical_transcriptome",main_id=self.option('main_table_id'))
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
        super(GenesetClassWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_run_log()
        self.set_db()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_geneset_go_class", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('medical_transcriptome.medical_transcriptome_geneset')
        self.logger.info("注释分类分析结果入库")
        if self.option("anno_type") == "go":
            time.sleep(30)
            output_file = self.option("geneset_go")

            ## 文件加入对应于基因id列的【基因名】这一新列
            add_name = GeneInfoSupple(task_id=self.option("task_id"), level=self.option("geneset_type"))
            try:
                geneid_list = [i for i in open(output_file).readline().strip().split('\t') if " list" in i]
                add_name.add_gene_name(output_file, ';', "go", geneid_list, output_file+'_addgenename')
                output_file__AddGenename = output_file+'_addgenename'
                print "geneid_list:    "+str(geneid_list)
            except:
                self.logger.info("无法增加【基因名】这一新列")
                print "无法增加【基因名】这一新列"
                output_file__AddGenename = output_file

            linkfile(output_file__AddGenename, self.output_dir + "/go_class_table.xls")
            api_geneset.add_go_regulate_detail(output_file, self.option("main_table_id"))
        print(output_file)
        self.end()

    def chart(self):
        has_chart = False
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"

        if self.option("anno_type") == "cog":
            cog_class_table = self.work_dir + "/cog_class_table.xls"
            count = len(open(cog_class_table, 'r').readlines())
            if count > 1:
                chart.chart_geneset_class_cog(cog_class_table)
                has_chart = True
        elif self.option("anno_type") == "go":
            go_class_table = self.work_dir + "/go_class_table.xls"
            count = len(open(go_class_table, 'r').readlines())
            if count>1:
                chart.chart_geneset_class_go(go_class_table)
                has_chart= True
        if has_chart:
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = glob.glob(self.work_dir + "/*.pdf")[0]
            linkfile(pdf_file, self.output_dir + "/{}_class.pdf".format(self.option("anno_type")))

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        linkfile(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("anno_type") == "go":
            self.inter_dirs = [
                ["04 GeneSet", "", "基因集分析结果目录",0],
                ["04 GeneSet/02 Annotation", "", "基因集功能注释", 0],
                ["04 GeneSet/02 Annotation/01 GO", "", "基因集GO功能注释", 0]
            ]
        if self.option("anno_type") == "go":
            result_dir.add_relpath_rules([
                [".", "", "基因集GO功能注释文件",0],
                ["go_class_table.xls", "xls", "基因集GO分类统计表 ",0],
                ['./run_parameter.txt', 'txt', '运行参数日志', 0],
                ["*.pdf", "pdf", "GO分类统计图",0],
            ])
        # print self.get_upload_files()
        super(GenesetClassWorkflow, self).end()
        # self.set_end()
        # self.fire('end')
        # self.end_unfinish_job()
        # self._upload_result()
        # self._import_report_data()
        # self._update("set_end")
        # self.step.finish()
        # self.step.update()
        # self.logger.info("运行结束!")
        # self._save_report_data()
        # super(GenesetClassWorkflow, self).end()

