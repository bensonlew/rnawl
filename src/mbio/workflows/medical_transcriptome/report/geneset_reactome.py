# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import glob
import shutil
import pandas as pd
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.ref_rna_v2.copy_file import CopyFile
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart_geneset import ChartGeneset
from mbio.packages.medical_transcriptome.gene_info_supple import GeneInfoSupple
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile


class GenesetReactomeWorkflow(Workflow):
    """
    基因集功能分类分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetReactomeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "reactome_annot", "type": "string"},
            {"name": "geneset_names", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "add_info", "type": "string"},
            {"name": "type", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "source", "type": "string"},
            {'name': 'reactome_version', 'type': 'string', 'default': "72"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.reactome_class = self.add_tool("medical_transcriptome.geneset.reactome_class")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/04 GeneSet/02 Annotation/03 Reactome')
        self.inter_dirs = []
        self.has_results = False
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete = InteractionDelete(bind_object=self, project_type="medical_transcriptome",
                                                  main_id=self.option('main_table_id'))
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
        super(GenesetReactomeWorkflow, self).send_log(data)

    def run(self):
        # super(GenesetReactomeWorkflow, self).run()
        self.reactome_class.on("end", self.set_db)
        self.get_run_log()
        self.run_reactome_class()
        super(GenesetReactomeWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_geneset_reactome_class", main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('medical_transcriptome.geneset_annot')
        if self.option("source") == "diff_exp":
            api_geneset = self.api.api('medical_transcriptome.diff_geneset_annot')
        self.logger.info("开始进行reactome_class的导表")
        pathway_file = self.reactome_class.output_dir + '/reactome_path.xls'
        class_file = self.reactome_class.output_dir + '/reactome_class.xls'
        svg_path = self.reactome_class.output_dir + '/svg'
        svg_gz_path = self.reactome_class.output_dir + '/svg.tar.gz'
        results_num = len(open(pathway_file, 'r').readlines())
        if results_num>1:
            self.has_results = True
            add_name = GeneInfoSupple(id2namedict = None, task_id = self.option("task_id"), level ="G")
            add_columns = ["Genes"]
            add_columns = [x + "_genes" for x in self.option("geneset_names").split(",")]

            try:
                add_name.add_gene_name(class_file, split=";", annot_type="reactome", add_columns=add_columns, outfile_path=self.output_dir + "/Reactome_class_stat.xls")
            except Exception as e:
                self.logger.info("添加基因name 失败 {}".format(e))
                CopyFile().linkfile(class_file, self.output_dir + "/Reactome_class_stat.xls")

            CopyFile().linkfile(pathway_file, self.output_dir + "/pathway_class_stat.xls")

            # CopyFile().linkdir(svg_path, self.output_dir + "/Reactome_pathways")
            CopyFile().linkfile(svg_gz_path, self.output_dir + "/Reactome_pathways.tar.gz")

            api_geneset.add_geneset_reactome_detail(self.option("main_table_id"),
                                            self.option("geneset_names"),
                                            pathway_file,
                                            class_file,
                                            svg_path,
                                            source=self.option('source'))


            # 更新主表
            conn = api_geneset.db["sg_geneset_reactome_class"]

            if self.option("source") == "diff_exp":
                conn = api_geneset.db["sg_diff_geneset_reactome_class"]
            self.workflow_output_tmp = self._sheet.output
            if re.match(r'tsanger:', self.workflow_output_tmp):
                self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
            else:
                self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
            graph_dir = os.path.join(self.workflow_output, 'svg')
            conn.update({"_id": ObjectId(self.option("main_table_id")) }, {"$set": {'graph_dir': graph_dir}}, upsert=True)
            # mark_file = glob.glob(self.output_dir + "/Reactome_pathways/*.mark")
            # for file in mark_file:
            #     os.remove(file)
        else:
            conn = api_geneset.db["sg_geneset_reactome_class"]
            conn.update({"_id": ObjectId(self.option("main_table_id"))}, {"$set": {'has_results': False}},
                        upsert=True)
        self.end()

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"
        if self.has_results:
            output_file = self.reactome_class.output_dir + '/reactome_class.xls'
            chart.chart_geneset_class_reactome(reactome_class_table=output_file)
            chart.to_pdf()
            # move pdf to result dir
            pdf_file = glob.glob(self.work_dir + "/*.pdf")
            for p in pdf_file:
                linkfile(p, self.output_dir + "/" + os.path.basename(p))
        else:
            pass

    # 这个没有去修改原来的tool,通过导表函数做了一些计算，所以上传文件来自work_dir
    def end(self):
        self.chart()
        ## set output
        '''
        mark_file = glob.glob(self.reactome_class.output_dir + "/pathways/*.html.mark")
        for file in mark_file:
            os.remove(file)
        '''
        # os.rename(self.reactome_class.output_dir + '/pathways', self.reactome_class.output_dir + '/reactome_pathways')
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        self.inter_dirs = [
            ["04 GeneSet", "", "基因集分析结果目录",0],
            ["04 GeneSet/02 Annotation", "", "基因集功能注释", 0],
            ["04 GeneSet/02 Annotation/03 Reactome", "", "基因集Reactome功能注释", 0]
        ]
        result_dir.add_relpath_rules([
            [".", "", "基因集Reactome功能注释文件",0],
            ["Reactome_class_stat.xls", "", "基因集Reactome分类统计表 ",0],
            ["pathway_class_stat.xls", "", "基因集Reactome通路统计表 ",0],
            ["Reactome_pathways.tar.gz", "", "基因集Reactome通路图",0],
            ["*.pdf", "pdf", "Reactome分类条形图", 0],
            ["Reactome_pathways", "", "基因集Reactome通路图",0],
            ['Reactome_pathways/*.svg', '', 'Reactome通路图片svg',0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0]
        ])
        result_dir.add_regexp_rules([
            [r'.*\.pdf', '', 'Reactome通路富集图片', 0],
        ])
        super(GenesetReactomeWorkflow, self).end()

    def run_reactome_class(self):
        opts = {
            'geneset_ids': self.option('geneset_id'),
            'reactome_annot': self.option('reactome_annot'),
            'reactome_version': self.option('reactome_version')
        }
        self.reactome_class.set_options(opts)
        self.reactome_class.run()
