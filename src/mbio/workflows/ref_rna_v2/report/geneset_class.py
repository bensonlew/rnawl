# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.workflow import Workflow
from collections import defaultdict
import os
import re
import time
from itertools import chain
from mbio.packages.ref_rna_v2.kegg_regulate import KeggRegulate
import json
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.interaction_rerun.interaction_delete import InteractionDelete,linkfile
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart_geneset import ChartGeneset
import glob


class GenesetClassWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetClassWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "geneset_go", "type": "string"},
            {"name": "geneset_cog", "type": "string"},
            {"name": "geneset_kegg", "type": "string"},
            {"name": "kegg_table", "type": "infile", "format": "ref_rna_v2.kegg_table"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "type", "type": "string"},
            {"name": "task_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "cog":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 GeneSet/02 COG_Annotation')
        if self.option("anno_type") == "go":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/03 GeneSet/03 GO_Annotation')

        self.inter_dirs = []
        if self._sheet.rerun:
            self.logger.info("该交互为重运行项目,先删除mongo库内容")
            interactiondelete =  InteractionDelete(bind_object=self,project_type="ref_rna_v2",main_id=self.option('main_table_id'))
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
        if self.option("anno_type") == "kegg":
            self.get_kegg_table()
        self.set_db()

    def get_run_log(self):
        if self.option("anno_type") == "cog":
            table = "sg_geneset_cog_class"
        elif self.option("anno_type") == "go":
            table = "sg_geneset_go_class"
        else:
            table = "sg_geneset_kegg_class"
        get_run_log = GetRunLog("ref_rna_v2", table=table, main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('ref_rna_v2.ref_rna_v2_geneset')
        self.logger.info("注释分类分析结果入库")
        if self.option("anno_type") == "cog":
            time.sleep(30)
            output_file = self.option("geneset_cog")
            # os.link(output_file, self.output_dir + "/cog_class_table.xls")
            linkfile(output_file,os.path.join(self.output_dir,"cog_class_table.xls"))
            api_geneset.add_geneset_cog_detail(output_file, self.option("main_table_id"))
        elif self.option("anno_type") == "go":
            time.sleep(30)
            output_file = self.option("geneset_go")
            # os.link(output_file, self.output_dir + "/go_class_table.xls")
            linkfile(output_file, os.path.join(self.output_dir, "go_class_table.xls"))
            api_geneset.add_go_regulate_detail(output_file, self.option("main_table_id"))
        else:
            output_file = self.output_dir + '/kegg_stat.xls'
            pathway_file = self.output_dir + '/pathways'
            api_geneset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
            api_geneset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        # os.link(output_file, self.output_dir + "/" + os.path.basename(output_file))
        print(output_file)
        self.end()

    def get_kegg_table(self):
        kegg = KeggRegulate()
        ko_genes, path_ko = self.option('kegg_table').get_pathway_koid()
        # regulate_dict = defaultdict(set)
        regulate_gene = {}
        with open(self.option("geneset_kegg"), "r") as f:
            for line in f:
                line = line.strip().split("\t")
                regulate_gene[line[0]] = line[1].split(",")  # multi_geneset_list
        pathways = self.output_dir + '/pathways'
        if not os.path.exists(pathways):
            os.mkdir(pathways)
        self.logger.info(ko_genes)
        self.get_kegg_pics(ko_genes=ko_genes, catgory=regulate_gene, out_dir=pathways)
        # kegg.get_pictrue(path_ko=path_ko, out_dir=pathways, regulate_dict=regulate_gene)  # 颜色
        # kegg.get_pictrue(path_ko=path_ko, out_dir=pathways)  # edited by shijin on 20170728
        kegg.get_regulate_table(ko_gene=ko_genes, path_ko=path_ko, regulate_gene=regulate_gene, output=self.output_dir + '/kegg_stat.xls')

    def get_kegg_pics(self, ko_genes, catgory, out_dir):
        kos = ko_genes.keys()
        genes = set()
        for ko in kos:
            genes.update(ko_genes[ko])
        gene_list = list(genes)
        color_dict = {}
        for gene in gene_list:
            color_dict[gene] = []
            # for catg in catgory.values():
            if len(catgory) == 1:
                if gene in catgory.values()[0]:
                    color_dict[gene].append("#00CD00")  # 绿色
            elif len(catgory) == 2:
                if gene in catgory.values()[0]:
                    color_dict[gene].append("#00CD00")  # 绿色
                elif gene in catgory.values()[1]:
                    color_dict[gene].append("#9932CC")  # 基佬紫
            else:
                pass
        ko_fg = {}
        for ko in kos:
            ko_fg[ko] = []
            for gene in ko_genes[ko]:
                ko_fg[ko]=ko_fg[ko]+color_dict[gene]
                ko_fg[ko] = list(set(ko_fg[ko]))
        with open(out_dir + "/kos", "w") as fw:
            fw.write("#KO\tbg\tfg\n")
            for ko in ko_fg.keys():
                fw.write(ko + "\tFF00FF\t" + ','.join(ko_fg[ko]) + "\n")
        self.logger.info("ko文件生成完毕")

        # pass

    def chart(self):
        chart = ChartGeneset()
        chart.work_dir = self.work_dir + "/"

        if self.option("anno_type") == "cog":
            cog_class_table = self.work_dir + "/cog_class_table.xls"
            chart.chart_geneset_class_cog(cog_class_table)
        elif self.option("anno_type") == "go":
            go_class_table = self.work_dir + "/go_class_table.xls"
            chart.chart_geneset_class_go(go_class_table)
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")[0]
        linkfile(pdf_file, self.output_dir + "/{}_class.pdf".format(self.option("anno_type")))
        # os.link(pdf_file, self.output_dir + "/{}_class.pdf".format(self.option("anno_type")))

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("anno_type") == "go":
            self.inter_dirs = [
                ["03 GeneSet", "", "基因集分析结果目录",0],
                ["03 GeneSet/03 GO_Annotation", "", "GO功能注释", 0]
            ]
        elif self.option("anno_type") == "cog":
            self.inter_dirs = [
                ["03 GeneSet", "", "基因集分析结果目录",0],
                ["03 GeneSet/02 COG_Annotation", "", "COG功能注释", 0]
            ]
        if self.option("anno_type") == "go":
            result_dir.add_relpath_rules([
                [".", "", "GO功能注释文件",0,"211526"],
                ["go_class_table.xls", "", "GO分类统计表",0,"211527"],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
                ["*.pdf", "pdf", "GO分类统计图",0],
            ])
        elif self.option("anno_type") == "cog":
            result_dir.add_relpath_rules([
                [".", "", "COG功能注释文件",0,"211528"],
                ["cog_class_table.xls", "", "COG分类统计表",0,"211529"],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
                ["*.pdf", "pdf", "COG分类统计图",0],
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

