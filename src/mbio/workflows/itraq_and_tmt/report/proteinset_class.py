# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last update liubinxu
from biocluster.workflow import Workflow
from collections import defaultdict
import os
import time
import re
from itertools import chain
from mbio.packages.denovo_rna.express.kegg_regulate import KeggRegulate
from mbio.packages.itraq_and_tmt.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import json
import glob


class ProteinsetClassWorkflow(Workflow):
    """
    蛋白集功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetClassWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "proteinset_go", "type": "string"},
            {"name": "proteinset_go2", "type": "string"},
            {"name": "proteinset_kegg", "type": "string"},
            {"name": "proteinset_cog", "type": "string"},
            {"name": "kegg_table", "type": "infile", "format": "denovo_rna_v2.kegg_table"},
            {"name": "anno_type", "type": "string"},
            {"name": "proteinset_type", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_table_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "proteinset_id", "type": "string"},
            {"name": "type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("anno_type") == "cog":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/03_Anno/03_COG')
        if self.option("anno_type") == "go":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/03_Anno/01_GO')
        if self.option("anno_type") == "go2":
            self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/5_Proteinset/03_Anno/01_GO')
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
        super(ProteinsetClassWorkflow, self).send_log(data)

    def run(self):
        self.start_listener()
        self.fire("start")
        if self.option("anno_type") == "kegg":
            self.get_kegg_table()
        self.set_db()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_proteinset = self.api.api('itraq_and_tmt.proteinset')
        self.logger.info("注释分类分析结果入库")
        if self.option("anno_type") == "cog":
            time.sleep(30)
            output_file = self.option("proteinset_cog")
            if os.path.exists(self.output_dir + "/cog_class_table.xls"):
                os.remove(self.output_dir + "/cog_class_table.xls")
            os.link(output_file, self.output_dir + "/cog_class_table.xls")
            api_proteinset.add_proteinset_cog_detail(output_file, self.option("main_table_id"))
        elif self.option("anno_type") == "go":
            time.sleep(30)
            output_file = self.option("proteinset_go")
            if os.path.exists(self.output_dir + "/go_class_table.xls"):
                os.remove(self.output_dir + "/go_class_table.xls")
            os.link(output_file, self.output_dir + "/go_class_table.xls")
            api_proteinset.add_go_regulate_detail(output_file, self.option("main_table_id"))
        elif self.option("anno_type") == "go2":
            output_file = self.option("proteinset_go2")
            if os.path.exists(self.output_dir + "/go_class_table.xls"):
                os.remove(self.output_dir + "/go_class_table.xls")
            os.link(output_file, self.output_dir + "/go_class_table.xls")
            api_proteinset.add_go_regulate_detail2(output_file, self.option("main_table_id"))
        else:
            output_file = self.output_dir + '/kegg_stat.xls'
            pathway_file = self.output_dir + '/pathways'
            api_proteinset.add_kegg_regulate_detail(self.option("main_table_id"), output_file)
            api_proteinset.add_kegg_regulate_pathway(pathway_file, self.option("main_table_id"))
        # os.link(output_file, self.output_dir + "/" + os.path.basename(output_file))
        print(output_file)
        self.end()

    def get_kegg_table(self):
        kegg = KeggRegulate()
        ko_genes, path_ko = self.option('kegg_table').get_pathway_koid()
        # regulate_dict = defaultdict(set)
        regulate_gene = {}
        with open(self.option("proteinset_kegg"), "r") as f:
            for line in f:
                line = line.strip().split("\t")
                regulate_gene[line[0]] = line[1].split(",")  # multi_proteinset_list
        pathways = self.output_dir + '/pathways'
        if not os.path.exists(pathways):
            os.mkdir(pathways)
        self.logger.info(ko_genes)
        self.get_kegg_pics(ko_genes=ko_genes, catgory=regulate_gene, out_dir=pathways)
        # kegg.get_pictrue(path_ko=path_ko, out_dir=pathways, regulate_dict=regulate_gene)  # 颜色
        # kegg.get_pictrue(path_ko=path_ko, out_dir=pathways)  # edited by shijin on 20170728
        # kegg.get_regulate_table(ko_gene=ko_genes, path_ko=path_ko, regulate_gene=regulate_gene, output=self.output_dir + '/kegg_stat.xls')

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
                if gene in catgory.values()[1]:
                    color_dict[gene].append("#00CD00")  # 绿色
            elif gene in catgory.values()[2]:
                    color_dict[gene].append("#9932CC")  # 基佬紫
            else:
                pass
        ko_fg = {}
        for ko in kos:
            ko_fg[ko] = []
            for gene in ko_genes[ko]:
                ko_fg[ko].expand(color_dict[gene])
                ko_fg[ko] = list(set(ko_fg[ko]))
        with open(out_dir + "/kos", "w") as fw:
            fw.write("#KO\tbg\tfg\n")
            for ko in ko_fg.keys():
                fw.write(ko + "\tFF00FF\t" + ko_fg[ko] + "\n")
        self.logger.info("ko文件生成完毕")

        # pass

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir+'/'
        chart.output_dir = self.work_dir+'/'
        if self.option("anno_type") in ["go","go2"]:
            proteinsetgo_class_file = os.path.join(self.output_dir, "go_class_table.xls")
            chart.chart_proteinsetcluster_goclass_web(proteinsetgo_class_file)
            chart.to_pdf()
            # move pdf to result dir
            for f in [['webmulti_go.go_bar.pdf','go_lev2.pdf'], ['webmulti_go_neg.bar_neg.pdf','go_lev2_both.pdf']]:
                if os.path.exists(os.path.join(self.work_dir, f[0])):
                    os.link(os.path.join(self.work_dir, f[0]), os.path.join(self.output_dir, f[1]))
        if self.option("anno_type") == "cog":
            cog_file = os.path.join(self.output_dir, "cog_class_table.xls")
            chart.chart_cog_web(cog_file)
            chart.to_pdf()
            # move pdf to result dir
            if os.path.exists(os.path.join(self.work_dir, "cog.bar.pdf")):
                os.link(os.path.join(self.work_dir, "cog.bar.pdf"), os.path.join(self.output_dir, "cog.pdf"))

    def end(self):
        self.chart()
        result_dir = self.add_upload_dir(self.output_dir)
        if self.option("anno_type") == "go":
            self.inter_dirs = [
                ["5_Proteinset", "", "蛋白集分析",0],
                ["5_Proteinset/03_Anno", "", "功能注释", 0],
                ["5_Proteinset/03_Anno/01_GO", "", "go注释", 0],
            ]
            result_dir.add_relpath_rules([
                [".", "", "蛋白集GO分类单向直方图结果目录", 0, "220072"],
                ["go_lev2.pdf", "", "GO分类统计柱形图", 0],
                ["go_lev2_both.pdf", "", "GO分类统计双向柱形图", 0],
                # ["./estimators.xls", "xls", "alpha多样性指数表"]
            ])
        if self.option("anno_type") == "go2":
            self.inter_dirs = [
                ["5_Proteinset", "", "蛋白集分析",0],
                ["5_Proteinset/03_Anno", "", "功能注释", 0],
                ["5_Proteinset/03_Anno/01_GO", "", "go注释", 0],
            ]
            result_dir.add_relpath_rules([
                [".", "", "蛋白集GO分类双向直方图结果目录", 0, "220073"],
                ["go_lev2.pdf", "", "GO分类统计柱形图", 0],
                ["go_lev2_both.pdf", "", "GO分类统计双向柱形图", 0],
                # ["./estimators.xls", "xls", "alpha多样性指数表"]
            ])
        elif self.option("anno_type") == "cog":
            self.inter_dirs = [
                ["5_Proteinset", "", "蛋白集分析",0],
                ["5_Proteinset/03_Anno", "", "功能注释", 0],
                ["5_Proteinset/03_Anno/03_COG", "", "cog注释", 0],
            ]
            result_dir.add_relpath_rules([
                [".", "", "蛋白集COG分类结果目录", 0, "220074"],
                ["cog.pdf", "", "COG注释统计图", 0],
            ])
        # print self.get_upload_files()
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
        super(ProteinsetClassWorkflow, self).end()
