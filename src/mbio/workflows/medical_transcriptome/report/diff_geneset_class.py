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
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class DiffGenesetClassWorkflow(Workflow):
    """
    基因集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffGenesetClassWorkflow, self).__init__(wsheet_object)
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
            self._sheet.output = self._sheet.output.replace('interaction_results',
                                                            'interaction_results/01 Diff_Express/03 DiffExpress_Geneset_Annotion/01 GO')

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
        super(DiffGenesetClassWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_run_log()
        if self.option("anno_type") == "kegg":
            self.get_kegg_table()
        self.set_db()

    def get_run_log(self):
        if self.option("anno_type") == "cog":
            table = "sg_diff_geneset_cog_class"
        elif self.option("anno_type") == "go":
            table = "sg_diff_geneset_go_class"
        else:
            table = "sg_diff_geneset_kegg_class"
        get_run_log = GetRunLog("medical_transcriptome", table=table, main_id=self.option('main_table_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()


    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_geneset = self.api.api('medical_transcriptome.diff_geneset')
        self.logger.info("注释分类分析结果入库")
        if self.option("anno_type") == "cog":
            time.sleep(30)
            output_file = self.option("geneset_cog")
            os.link(output_file, self.output_dir + "/cog_class_table.xls")
            api_geneset.add_geneset_cog_detail(output_file, self.option("main_table_id"))
        elif self.option("anno_type") == "go":
            time.sleep(30)
            output_file = self.option("geneset_go")
            os.link(output_file, self.output_dir + "/go_class_table.xls")
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
        # self.logger.info(ko_genes)
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

    def end(self):
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("anno_type") == "go":
            self.inter_dirs = [
                ["01 Diff_Express", "", "差异基因数据挖掘结果目录", 0],
                ["01 Diff_Express/03 DiffExpress_Geneset_Annotion", "", "差异基因集功能注释分析", 0],
                ["01 Diff_Express/03 DiffExpress_Geneset_Annotion/01 GO", "", "差异基因集GO功能注释分析", 0]
            ]
            result_dir.add_relpath_rules([
                [".", "", "差异基因集DO功能注释文件",0],
                ["go_class_table.xls", "", "GO分类统计表",0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])

        # print self.get_upload_files()
        super(DiffGenesetClassWorkflow, self).end()
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
        # super(DiffGenesetClassWorkflow, self).end()

