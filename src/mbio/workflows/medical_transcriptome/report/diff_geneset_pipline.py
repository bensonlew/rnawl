# -*- coding: utf-8 -*-
# __author__ = 'fwy'

from biocluster.workflow import Workflow
import types
from bson.objectid import ObjectId
import os
import re
import glob
import shutil
import pandas as pd
import gevent
from collections import OrderedDict
import unittest
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
# from mbio.packages.medical_transcriptome.upload import set_diff_pipline
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.upload import Upload
from mbio.packages.medical_transcriptome.copy_file import CopyFile



class DiffGenesetPiplineWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffGenesetPiplineWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "diff_id", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "species", "type": "string"},
            {"name": "pipline_main_id", "type": "string"},
            {"name": "level", "type": "string"},
            {"name": "genesets", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "submit_location", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Diff_Express')
        self.inter_dirs = []
        self.genset_dict= OrderedDict()
        self.file_prepare = self.add_module("medical_transcriptome.diff_geneset.diff_geneset_prepare")
        self.chart = self.add_tool("medical_transcriptome.chart")
        self.geneset_analysis = []
        self.kegg_level_path = ""
        self.geneset_list=""
        self.all_list = ""
        self.analysis_names =[]

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
        super(DiffGenesetPiplineWorkflow, self).send_log(data)

    def run(self):

        self.get_run_log()
        self.run_file_prepare()
        super(DiffGenesetPiplineWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_diff_geneset_pipline", main_id=self.option('pipline_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_file_prepare(self):
        opts = {
            "geneset_names": self.option("genesets"),
            "diff_id" : self.option("diff_id"),
            "level" : self.option("level"),
            "task_id" :self.option("task_id")
        }
        self.file_prepare.set_options(opts)
        self.file_prepare.on("end",self.run_genesets_analysis)
        self.file_prepare.run()

    def run_genesets_analysis(self):
        file_json_path = os.path.join(self.file_prepare.output_dir,"prepare_json")
        if not os.path.exists(file_json_path):
            raise Exception("未找到文件准备汇总json文件")
        with open(file_json_path,"r") as j:
            file_dict = json.load(j)
        self.kegg_level_path = file_dict["common_file"]["common_annot_file"]["kegg_level_table"]

        gensets_infos = file_dict["genesets"]
        for geneset_info in sorted(gensets_infos.keys()):
            opts = self.opts_prepare(file_dict,gensets_infos[geneset_info])
            geneset_analysis = self.add_module("medical_transcriptome.diff_geneset.diff_geneset_analysis")
            geneset_analysis.set_options(opts)
            self.geneset_analysis.append(geneset_analysis)
        if self.geneset_analysis:
            if len(self.geneset_analysis) > 1:
                self.on_rely(self.geneset_analysis, self.set_output)
            elif len(self.geneset_analysis) == 1:
                self.geneset_analysis[0].on('end', self.set_output)
        else:
            self.set_error("geneset_analysis列表为空！")
        for tool in self.geneset_analysis:
            gevent.sleep(1)
            tool.run()


    def opts_prepare(self,file_dict,genesets):
        self.analysis_names = ["go","kegg"]
        if self.option("level") == "G":
            self.analysis_names.append("reactome")
            if self.option("species") == "Homo_sapiens":
                self.analysis_names.append("do")
        if genesets["gene_num"] == 0:
            opts = {
                "task_id": self.option("task_id"),
                "level": self.option("level"),
                "species": self.option("species"),
                "genes_num": genesets["gene_num"],
                "geneset_name": genesets["geneset_name"],
                "geneset_path": genesets["geneset_path"],
                "regulate": genesets["regulate"],
            }
        else :
            opts={
                "task_id":self.option("task_id"),
                "level" :self.option("level"),
                "species": self.option("species"),
                "go_list":file_dict["common_file"]["common_annot_file"]["go_list"],
                "kegg_table":file_dict["common_file"]["common_annot_file"]["kegg_table"],
                "kegg_table2": file_dict["common_file"]["common_annot_file"]["kegg_level_table"],
                "all_list": file_dict["common_file"]["common_annot_file"]["all_list"],
                "add_info" : file_dict["common_file"]["common_annot_file"]["add_info"],
                "kegg_version" : file_dict["common_file"]["common_annot_file"]["kegg_version"],
                "geneset_name": genesets["geneset_name"],
                "genes_num": genesets["gene_num"],
                "geneset_path":genesets["geneset_path"],
                "regulate": genesets["regulate"],
                "gene_list":genesets["file_path"]["go_enrich"]["gene_list_path"],
                "gene_multi" : genesets["file_path"]["kegg_class"]["multi_gene_list_path"],
                "go_class" :genesets["file_path"]["go_class"]["go_class_path"]
            }
            if "reactome" in self.analysis_names:
                opts.update({
                    "reactome_annot":file_dict["common_file"]["common_annot_file"]["reactome_annot"],
                    "reactome_version": file_dict["common_file"]["common_annot_file"]["reactome_version"],
                })
            if "do" in self.analysis_names:
                opts.update({
                    "do_list":file_dict["common_file"]["common_annot_file"]["do_list"],
                    "do_class": genesets["file_path"]["do_class"]["do_class_path"],
                })
        return opts

    def set_output(self):
        for geneset_analysis in self.geneset_analysis:
            analysis_json = os.path.join(geneset_analysis.output_dir,"analysis_json")
            with open(analysis_json,"r") as f:
                analysis_dict = json.load(f)
            geneset_name = analysis_dict["geneset_name"]
            genset_result_dir = os.path.join(self.output_dir,geneset_name)
            if os.path.exists(genset_result_dir):
                shutil.rmtree(genset_result_dir)

            # add by fwy 20210125
            CopyFile().linkdir(geneset_analysis.output_dir, genset_result_dir)
            # os.system('ln -s {} {}'.format(geneset_analysis.output_dir,genset_result_dir))
            # shutil.copytree(geneset_analysis.output_dir,genset_result_dir)
        self.logger.info("准备导表喽")
        self.run_chart()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        self.logger.info("开始导表喽")
        if os.path.exists(os.path.join(self.work_dir,"temporary")):
            shutil.rmtree(os.path.join(self.work_dir,"temporary"))
        os.makedirs(os.path.join(self.work_dir,"temporary"))
        self.export_temporary = os.path.join(self.work_dir,"temporary")
        api_geneset_pipline = self.api.api('medical_transcriptome.diff_geneset_pipline')
        api_geneset_pipline.add_diff_genest_pipline_table(self.output_dir,diff_geneset_pipline_id = self.option("pipline_main_id"),
                                                          diff_id=self.option("diff_id"),task_id =self.option("task_id"),
                                                          analysis_names = self.analysis_names,kegg_level_path =self.kegg_level_path,inter_path=self.export_temporary
                                                          )
        self.logger.info("导表导完喽")
        self.modify_upload()

    def modify_upload(self):
        self.logger.info("调整文件格式喽")
        upload_dir = os.path.join(self.work_dir, "upload")
        upload = Upload(task_id=self.option("task_id"),level =self.option("level"))
        upload.set_diff_pipline(idir=self.output_dir, odir=upload_dir, species=self.option("species"),
                                level=self.option("level"))
        self.logger.info("文件格式调整喽")
        self.logger.info("准备生成图片喽")
        self.end()

    def run_chart(self):
        self.logger.info("开始生成图片喽")
        genesets = os.listdir(self.output_dir)
        chart_dict = {
            "type": "diff_pipline",
            "genesets": genesets,
            "go_class":"{table_dir}/{geneset_name}/diff_go_class/go_class_table.xls".format(table_dir = self.output_dir,geneset_name = "{geneset_name}"),
            "go_enrich": "{table_dir}/{geneset_name}/diff_go_enrich/go_enrich_geneset_list_gene.xls".format(
                table_dir=self.output_dir, geneset_name="{geneset_name}"),
            "kegg_level": self.kegg_level_path,
            "kegg_class": "{table_dir}/{geneset_name}/diff_kegg_class/kegg_stat.xls".format(
                table_dir=self.output_dir, geneset_name="{geneset_name}"),
            "kegg_enrich": "{table_dir}/{geneset_name}/diff_kegg_enrich/enrich/{geneset_name}_gene.list.DE.list.check.kegg_enrichment.xls".format(
                table_dir=self.output_dir, geneset_name="{geneset_name}")
        }
        if "reactome" in self.analysis_names:
            chart_dict.update({
                "reactome_class": "{table_dir}/{geneset_name}/diff_reactome_class/reactome_class.xls".format(
                    table_dir=self.output_dir, geneset_name="{geneset_name}"),
                "reactome_enrich": "{table_dir}/{geneset_name}/diff_reactome_enrich/{geneset_name}_gene.list.reactome_enrichment.xls".format(
                    table_dir=self.output_dir, geneset_name="{geneset_name}")
            })
        if "do" in self.analysis_names:
            chart_dict.update({
                "do_class": "{table_dir}/{geneset_name}/diff_do_class/do_level2_class.xls".format(
                    table_dir=self.output_dir, geneset_name="{geneset_name}"),
                "do_enrich": "{table_dir}/{geneset_name}/diff_do_enrich/do_enrichment.xls".format(
                    table_dir=self.output_dir, geneset_name="{geneset_name}")
            })


        with open(self.work_dir + "/chart_workflow.json", 'w') as json_f:
            json.dump(chart_dict, json_f, sort_keys=True, indent=4)
        self.chart.set_options({
            "file_json": self.work_dir + "/chart_workflow.json"
        })
        self.chart.on('end', self.set_db)
        self.chart.run()

        # if "reactome" in self.analysis_names:
        #     chart_dict.update({
        #         "reactome_class": "{table_dir}/{geneset_name}/diff_go_class/go_class_table.xls".format(
        #             table_dir=self.output_dir, geneset_name="{geneset_name}"),
        #         "reactome_enrich": "{table_dir}/{geneset_name}/diff_go_class/go_class_table.xls".format(
        #             table_dir=self.output_dir, geneset_name="{geneset_name}"),
        #         "go_class": "{output_dir}/{geneset_name}/diff_go_class/go_class_table.xls".format(
        #             output_dir=self.output_dir, geneset_name="{geneset_name}"),
        #         "go_class": "{output_dir}/{geneset_name}/diff_go_class/go_class_table.xls".format(
        #             output_dir=self.output_dir, geneset_name="{geneset_name}"),
        #     })

    def move_chart_file(self):
        upload_dir = os.path.join(self.work_dir, "upload")
        genesets = os.listdir(self.output_dir)
        file2uploads = [
            ("{}.go_annot.gene_set.column.pdf".format("{geneset_name}"),
             "03DiffExpress_Geneset_Annotion/01GO/{}".format("{geneset_name}")),
            ("{}.go_enrich.gene_set.*.pdf".format("{geneset_name}"),
             "04DiffExpress_Geneset_Enrich/01GO/{}".format("{geneset_name}")),
            ("{}*kegg_annot.gene_set.column.pdf".format("{geneset_name}"),
             "03DiffExpress_Geneset_Annotion/02KEGG/{}".format("{geneset_name}")),
            ("{}.kegg_enrich.gene_set.*.pdf".format("{geneset_name}"),
             "04DiffExpress_Geneset_Enrich/02KEGG/{}".format("{geneset_name}")),
            ("{}*.reactome_annot.gene_set.column.pdf".format("{geneset_name}"),
             "03DiffExpress_Geneset_Annotion/03Reactome/{}".format("{geneset_name}")),
            ("{}.reactome_enrich.gene_set.*.pdf".format("{geneset_name}"),
             "04DiffExpress_Geneset_Enrich/03Reactome/{}".format("{geneset_name}")),
            ("{}.do_annot.gene_set.column.pdf".format("{geneset_name}"),
             "03DiffExpress_Geneset_Annotion/04DO/{}".format("{geneset_name}")),
            ("{}.do_enrich.gene_set.*.pdf".format("{geneset_name}"),
             "04DiffExpress_Geneset_Enrich/04DO/{}".format("{geneset_name}"))
        ]
        for geneset in genesets :
            for filefrom, fileto in file2uploads:
                pdf_file = glob.glob((self.chart.work_dir + "/" + filefrom).format(geneset_name = geneset))
                for p in pdf_file:
                    if os.path.exists(upload_dir + "/" + fileto.format(geneset_name = geneset) + "/" + os.path.basename(p)):
                        os.remove(upload_dir + "/" +  fileto.format(geneset_name = geneset) + "/" + os.path.basename(p))
                    os.link(p, upload_dir + "/" +  fileto.format(geneset_name = geneset) + "/" + os.path.basename(p))
        pass


    def end(self):
        self.move_chart_file()
        upload_dir = os.path.join(self.work_dir, "upload")
        result_dir = self.add_upload_dir(upload_dir)
        self.inter_dirs = [
            ["01 Diff_Express", "", "差异基因数据挖掘结果目录", 0],
        ]
        result_dir.add_relpath_rules([
            [r'.', '', '差异基因一键化结果目录', 0],
            ["03DiffExpress_Geneset_Annotion", "", "差异基因集功能注释分析", 0],
            ["04DiffExpress_Geneset_Enrich", "", "差异基因集功能富集分析", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        result_dir.add_regexp_rules([
            [r'03DiffExpress_Geneset_Annotion/01GO/.*', '', '差异基因集GO功能注释分析', 0],
            [r'03DiffExpress_Geneset_Annotion/01GO/.*/go_class_stat\.xls', '', 'GO分类统计表 ', 0],
            [r'03DiffExpress_Geneset_Annotion/01GO/.*/go_detail\.xls ', '', '基因/转录本对应GO注释详情表', 0],
            [r'03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_class_stat\.xls', '', 'KEGG分类统计表', 0],
            [r'03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_pathways', '', 'KEGG通路图', 0],
            [r'03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_pathways/.*png', '', 'KEGG通路图片png', 0],
            [r'03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_pathways/.*html', '', 'KEGG通路html文件', 0],
            [r'03DiffExpress_Geneset_Annotion/03Reactome/.*', '', '差异基因集功能注释分析', 0],
            [r'03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_class_stat\.xls', '',
             'Reactome分类统计表 ', 0],
            [r'03DiffExpress_Geneset_Annotion/03Reactome/.*/pathway_class_stat\.xls', '',
             'Reactome通路统计表 ', 0],
            [r'03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_pathways', '', 'Reactome通路图 ', 0],
            [r'03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_pathways/.*\.png', '',
             'Reactome通路图片png', 0],
            [r'03DiffExpress_Geneset_Annotion/03Reactome/.*/Reactome_pathways/.*\.html', '',
             'Reactome通路html文件 ', 0],
            [r'03DiffExpress_Geneset_Annotion/04DO/.*', '', '差异基因集功能注释分析', 0],
            [r'03DiffExpress_Geneset_Annotion/04DO/.*/DO_class_stat\.xls ', '', 'DO分类统计表 ', 0],
            [r'04DiffExpress_Geneset_Enrich/01GO/.*', '', '差异基因集功能富集分析', 0],
            [r'04DiffExpress_Geneset_Enrich/01GO/.*/go_enrich_stat\.xls', '', 'GO富集分析统计表  ', 0],
            [r'04DiffExpress_Geneset_Enrich/01GO/.*/go_detail\.xls ', '', '基因/转录本对应GO富集详情表', 0],
            [r'04DiffExpress_Geneset_Enrich/02KEGG/.*/kegg_enrich_stat\.xls', '', 'KEGG富集分析统计表 ', 0],
            [r'04DiffExpress_Geneset_Enrich/02KEGG/.*', '', '差异基因集功能富集分析 ', 0],
            [r'04DiffExpress_Geneset_Enrich/02KEGG/.*/kegg_enrich_pathways', '', 'KEGG富集通路图', 0],
            [r'04DiffExpress_Geneset_Enrich/02KEGG/.*/kegg_enrich_pathways/.*png', '', 'KEGG通路图片png', 0],
            [r'04DiffExpress_Geneset_Enrich/02KEGG/.*/kegg_enrich_pathways/.*html', '', 'KEGG通路html文件',
             0],
            [r'04DiffExpress_Geneset_Enrich/03Reactome/.*', '', '差异基因集功能富集分析', 0],
            [r'04DiffExpress_Geneset_Enrich/03Reactome/.*/Reactome_enrich_stat\.xls', '',
             ' Reactome富集分析统计表  ', 0],
            [r'04DiffExpress_Geneset_Enrich/03Reactome/.*/Reactome_pathways', '', 'Reactome通路图 ', 0],
            [r'04DiffExpress_Geneset_Enrich/03Reactome/.*/Reactome_pathways/.*\.png', '',
             'Reactome通路图片png', 0],
            [r'04DiffExpress_Geneset_Enrich/03Reactome/.*/Reactome_pathways/.*\.html', '',
             'Reactome通路html文件 ', 0],
            [r'04DiffExpress_Geneset_Enrich/04DO/.*', '', '差异基因集功能富集分析', 0],
            [r'04DiffExpress_Geneset_Enrich/04DO/.*/DO_enrich_stat\.xls ', '', 'DO富集分析统计表 ', 0],
            #新增图片的文件描述
            [r"03DiffExpress_Geneset_Annotion/01GO/.*/.*go_annot\.gene_set\.column\.pdf", 'pdf',
             'GO注释柱状图', 0],
            [r"04DiffExpress_Geneset_Enrich/01GO/.*/.*bar\.pdf", 'pdf', 'GO富集图片(柱形图)', 0],
            [r"04DiffExpress_Geneset_Enrich/01GO/.*/.*bar_line\.pdf", 'pdf', 'GO富集图片(柱形图-带折线)', 0],
            [r"04DiffExpress_Geneset_Enrich/01GO/.*/.*buble\.pdf", 'pdf', 'GO富集图片(气泡图)', 0],
            [r"04DiffExpress_Geneset_Enrich/01GO/.*/.*buble2\.pdf", 'pdf', 'GO富集图片(气泡图-分散型)', 0],
            [r"03DiffExpress_Geneset_Annotion/02KEGG/.*/kegg_annot\.{}*genes\.column\.pdf", 'pdf',
             'kegg注释柱状图', 0],
            [r"04DiffExpress_Geneset_Enrich/02KEGG/.*/.*bar\.pdf", 'pdf', 'kegg富集相关图片(柱形图)', 0],
            [r"04DiffExpress_Geneset_Enrich/02KEGG/.*/.*bar_line\.pdf", 'pdf', 'kegg富集相关图片(柱形图-带折线)', 0],
            [r"04DiffExpress_Geneset_Enrich/02KEGG/.*/.*buble\.pdf", 'pdf', 'kegg富集相关图片(气泡图)', 0],
            [r"04DiffExpress_Geneset_Enrich/02KEGG/.*/.*buble2\.pdf", 'pdf', 'kegg富集相关图片(气泡图-分散型)', 0],
            [r"03DiffExpress_Geneset_Annotion/03Reactome/.*/.*reactome_annot\.gene_set\.column.pdf",
             'pdf',
             'reactome注释柱状图', 0],
            # [r"04DiffExpress_Geneset_Enrich/03Reactome/.*/.*reactome_enrich\.gene_set.*\.pdf", 'pdf',
            #  'reactome富集相关图片', 0],
            [r"04DiffExpress_Geneset_Enrich/03Reactome/.*/.*bar\.pdf", 'pdf', 'reactome富集相关图片(柱形图)', 0],
            [r"04DiffExpress_Geneset_Enrich/03Reactome/.*/.*bar_line\.pdf", 'pdf',
             'reactome富集相关图片(柱形图-带折线)', 0],
            [r"04DiffExpress_Geneset_Enrich/03Reactome/.*/.*buble\.pdf", 'pdf', 'reactome富集相关图片(气泡图)',
             0],
            [r"04DiffExpress_Geneset_Enrich/03Reactome/.*/.*buble2\.pdf", 'pdf',
             'reactome富集相关图片(气泡图-分散型)', 0],
            [r"03DiffExpress_Geneset_Annotion/04DO/.*/.*do_annot\.gene_set\.column\.pdf", 'pdf',
             'DO注释柱状图', 0],
            # [r"04DiffExpress_Geneset_Enrich/04DO/.*/.*do_enrich\.gene_set.*\.pdf", 'pdf',
            #  'DO富集相关图片', 0],
            [r"04DiffExpress_Geneset_Enrich/04DO/.*/.*bar\.pdf", 'pdf', 'DO富集相关图片(柱形图)', 0],
            [r"04DiffExpress_Geneset_Enrich/04DO/.*/.*bar_line\.pdf", 'pdf', 'DO富集相关图片(柱形图-带折线)', 0],
            [r"04DiffExpress_Geneset_Enrich/04DO/.*/.*buble\.pdf", 'pdf', 'DO富集相关图片(气泡图)', 0],
            [r"04DiffExpress_Geneset_Enrich/04DO/.*/.*buble2\.pdf", 'pdf', 'DO富集相关图片(气泡图-分散型)', 0],
            [r"01Diff_Express/05DiffExp_Venn/diff_genesets\.analysis\.venn\.pdf", 'pdf',
             '差异基因集venn分析', 0],
            [r"01Diff_Express/02DiffExpress_Cluster_Analysis/.*cluster.*\.pdf", 'pdf',
             '差异基因集聚类分析图片', 0]
            
        ])
        super(DiffGenesetPiplineWorkflow, self).end()

