# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.config import Config
import re
import collections
import json
from itertools import islice
import subprocess
import gridfs
import os
import sys
import pandas as pd
import json
import math
import numpy as np
from scipy import stats
import xml.etree.ElementTree as ET
import lxml.html
from collections import OrderedDict
from chart import Chart
from bson.objectid import ObjectId
import glob

class ChartGeneset(Chart):
    def __init__(self):
        super(ChartGeneset, self).__init__()


    def chart_geneset_venn(self, geneset_ids):
        project_type = 'denovo_rna_v2'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        collection = db['sg_geneset_detail']
        main_collection = db['sg_geneset']

        source = list()
        for geneset_id in geneset_ids.split(","):
            my_result = main_collection.find_one({'main_id': ObjectId(geneset_id)})
            name = my_result["name"]
            print name
            results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
            seq_list = results["seq_list"]

            source.append({
                "data": seq_list,
                "name": name
            })

        self.chart_venn("geneset", "",  source, "dgeneset.venn.venn.json")
        # if len(source) <= 6:
        #     self.chart_upset("geneset", "",  source, "geneset.venn.upset.json")

    def chart_geneset_cluster(self, cluster_exp, cluster_tree, sample_tree, subcluster_list, group_dict=None, samples_order=None):
        cluster_pd = pd.read_table(cluster_exp, index_col=0, header=0)
        cluster_pd.index = cluster_pd.index.map(lambda x:x.replace("(", "-").replace(")", "-").replace(":", "-"))
        # samples 聚类图样本顺序， samples_order 折线图样本顺序
        if os.path.exists(sample_tree):
            with open(sample_tree, 'r') as f:
                sample_tree = f.readline().strip()
                samples = f.readline().strip().split(";")
                # samples_order = samples
        else:
            sample_tree = None
            samples = list(cluster_pd.columns)

        if os.path.exists(cluster_tree):
            with open(cluster_tree, 'r') as f:
                gene_tree = f.readline().strip()
                genes = f.readline().strip().split(";")
        else:
            gene_tree = None
            genes = list(cluster_pd.index)


        # 折线按表达量表排序
        if not samples_order:
            samples_order = list(cluster_pd.columns)

        cluster_pd["seq_id"] = cluster_pd.index
        cluster_pd = cluster_pd.loc[genes, :]
        corr_heat = [[""] + samples]
        for line_dict in cluster_pd.iterrows():
            corr_heat.append(
                [line_dict[0]] + [line_dict[1][s] for s in samples]
            )

        sample2group_source = [["name", "group"]]
        if group_dict:
            sample2group_dict = dict()
            for g,ss in group_dict.items():
                for s in ss:
                    # sample2group_dict
                    sample2group_source.append([s, g])
        else:
            for s in samples:
                sample2group_source.append([s, s])

        gene2group_source = [["name", "group"]]
        gene2group_dict = dict()
        for c in subcluster_list:
            print c
            c_num = c.split("subcluster_")[1].split(".")[0].split("_")[0]
            sub_genes = list()
            with open(c, 'r') as f:
                f.readline()
                for line in f:
                    gid = line.strip().split()[0]
                    gid = gid.replace("(", "-").replace(")", "-").replace(":", "-")
                    sub_genes.append(gid)
                    gene2group_dict[gid] = int(c_num)
            # print c

            source_line = list()
            sub_pd = cluster_pd.loc[sub_genes, :]
            for rec in sub_pd.to_dict("records"):
                # print rec
                source_line.append({
                    "data": [rec[s] for s in samples_order],
                    "group": "protein",
                    "name": rec["seq_id"]
                })
            source_line.append({
                "data": [sub_pd[s].mean() for s in samples_order],
                "group": "mean",
                "name": "mean"
            })
            # print source_line[:3]
            source_point = [["name", "x", "y", "category"]]
            for s in samples:
                source_point.append(["mean", s, sub_pd[s].mean(), "mean"])

            title = "subcluster_{}({} genes)".format(c_num, len(sub_pd))
            categories = samples_order
            self.chart_line_point("cluster", "." + c_num, source_line, source_point, categories, title, "geneset.cluster.line.json")

        for gene in genes:
            if gene in gene2group_dict:
                gene2group_source.append([gene, gene2group_dict[gene]])

        self.chart_heat_tree("geneset", ".cluster", corr_heat, sample_tree, gene_tree,  sample2group_source, gene2group_source, "geneset.cluster.heatmap.json")

    def chart_geneset_gsva_cluster(self, cluster_exp, cluster_tree, sample_tree, group_dict=None, samples_order=None):
        with open(sample_tree, 'r') as f:
            sample_tree = f.readline().strip()
            samples = f.readline().strip().split(";")

        with open(cluster_tree, 'r') as f:
            gene_tree = f.readline().strip()
            genes = f.readline().strip().split(";")

        cluster_pd = pd.read_table(cluster_exp, index_col=0, header=0)

        cluster_pd["seq_id"] = cluster_pd.index
        cluster_pd = cluster_pd.loc[genes, :]
        corr_heat = [[""] + samples]
        for line_dict in cluster_pd.iterrows():
            corr_heat.append(
                [line_dict[0]] + [line_dict[1][s] for s in samples]
            )

        sample2group_source = [["name", "group"]]
        if group_dict:
            for g,ss in group_dict.items():
                for s in ss:
                    sample2group_source.append([s, g])
        else:
            for s in samples:
                sample2group_source.append([s, s])

        gene2group_source = [["name", "group"]]
        for gene in genes:
            gene2group_source.append([gene, gene])

        self.chart_heat_tree("geneset_gsva", ".cluster", corr_heat, sample_tree, gene_tree,  sample2group_source, gene2group_source, "geneset.gsva.cluster.heatmap.json")

    def chart_geneset_gsva_diff(self, diff_file):
        name = os.path.basename(diff_file).split('.limma')[0]
        limma_file = pd.read_table(diff_file, header=0, index_col=0, sep='\t')
        limma_file.sort_values(by='padjust', inplace=True)
        if limma_file.shape[0] < 20:
            limma_20 = limma_file
        else:
            limma_20 = limma_file[0:20]
        limma_20.sort_values(by='log2fc', ascending=False, inplace=True)
        geneset_list = limma_20.index.tolist()
        visualColorValue = list()
        visualColorValue_list = {'up': "#FF0000", 'nosig': "#808080", 'down': "#008000"}
        source = list()
        source.append(["item", "value", "type"])
        for i in geneset_list:
            if limma_20.loc[i]['padjust'] >= 0.05:
                significant = 'nosig'
                source.append([i, limma_20.loc[i]['log2fc'], significant])
            else:
                if limma_20.loc[i]['log2fc'] >= 0.2:
                    significant = 'up'
                    source.append([i, limma_20.loc[i]['log2fc'], significant])
                if -0.2 < limma_20.loc[i]['log2fc'] < 0.2:
                    significant = 'nosig'
                    source.append([i, limma_20.loc[i]['log2fc'], significant])
                if limma_20.loc[i]['log2fc'] <= -0.2:
                    significant = 'down'
                    source.append([i, limma_20.loc[i]['log2fc'], significant])
            if visualColorValue_list[significant] not in visualColorValue:
                visualColorValue.append(visualColorValue_list[significant])

        # up_dict = {'data': [], 'name': 'up'}
        # nosig_dict = {'data': [], 'name': 'nosig'}
        # down_dict = {'data': [], 'name': 'down'}
        # for i in geneset_list:
        #     if limma_20.loc[i]['padjust'] >= 0.05:
        #         nosig_dict['data'].append(limma_20.loc[i]['log2fc'])
        #     else:
        #         if limma_20.loc[i]['log2fc'] >= 0.2:
        #             up_dict['data'].append(limma_20.loc[i]['log2fc'])
        #         if -0.2 < limma_20.loc[i]['log2fc'] < 0.2:
        #             nosig_dict['data'].append(limma_20.loc[i]['log2fc'])
        #         if limma_20.loc[i]['log2fc'] <= -0.2:
        #             down_dict['data'].append(limma_20.loc[i]['log2fc'])
        self.chart_gsva_diff('geneset_gsva_diff.', name, source, visualColorValue, 'gsva.diff.json')

    def chart_geneset_class_cog(self, cog_class_table, geneset_list=None):
        a = pd.read_table(cog_class_table, header=0, index_col=0)
        col_names = a.columns
        b = a.iloc[:, :-1]
        b.columns = col_names[1:]
        b.sort_values(by='Functional Categoris', ascending=True, inplace=True)
        categories = [x.split()[0][1] for x in list(b["Functional Categoris"])]
        if geneset_list:
            pass
        else:
            geneset_list = [c.split("_COG")[0]  for c in b.columns if c.endswith("_COG")]

        source = list()
        for geneset in geneset_list:
            source.append({
                "data": list(b[geneset + "_COG"]),
                "name": geneset
            })
        self.chart_class_column("cog_annot", ".gene_set", source, categories, None, None, "geneset.annot_cog_stat.bar.json")

    def chart_geneset_class_go(self, go_class_table, geneset_list=None, top=20):
        a = pd.read_table(go_class_table, header=0, index_col=0)
        # 原表格多一列处理
        col_names = a.columns
        b = a.iloc[:, :-1]
        b.columns = col_names[1:]
        b["sum"] = sum([b[x] for x in col_names if x.endswith("num")])
        c = b.sort_values(by=["sum"], ascending=[False])[:top]
        d = c.sort_index(ascending=False)
        d["type"] = d.index
        d = d.sort_values(by=["type", "sum"], ascending=[False, False])
        categories = list(d["Term"])
        if geneset_list:
            pass
        else:
            geneset_list = [c.split(" ")[0] for c in a.columns if c.endswith("num")]

        source = list()

        if len(geneset_list) == 1:
            # 单基因集
            geneset = geneset_list[0]
            source = [
                ["item"] + list(d["Term"]),
                ["series"] + list(d[geneset + " num"]),
                ["category"] + list(d["type"])
            ]
        else:
            for geneset in geneset_list:
                try:
                    percent = [float(p.split("(")[0]) for p in d[geneset + " percent"]]
                except:
                    percent = [float(0) for p in d[geneset + " percent"]]
                source.append({
                    "percent": percent,
                    "data": list(d[geneset + " num"]),
                    "name": geneset
                })


        class_source = [["p1", "p2", "categoy"]]
        for atype,sub_class in d.groupby("type", sort=False):
            class_source.append([sub_class.iloc[0,]["Term"], sub_class.iloc[-1,]["Term"], atype])
        if len(geneset_list) >= 2:
            bac_json = "geneset.annot_go_stat.bar.json"
        else:
            bac_json = "geneset.annot_go_stat2.bar.json"

        self.chart_class_column("go_annot", ".gene_set", source, categories, class_source, None, bac_json)

    def chart_geneset_class_kegg(self, kegg_class_table, geneset_list=None, top=20):
        a = pd.read_table(kegg_class_table, header=0)
        # a = a[:top]

        geneset_list = [c.split("_num")[0] for c in a.columns if c.endswith("num")]
        for geneset in geneset_list:
            a_choose = a[a[geneset + "_num"]>0]
            source = [
                ["item"] + list(a_choose["second_category"]),
                ["series"] + list(a_choose[geneset + "_num"]),
                ["category"] + list(a_choose["first_category"])
            ]

            class_source = ["p1", "p2", "category"]
            for atype, sub_class in a.groupby("first_category"):
                class_source.append([sub_class.iloc[0,]["second_category"], sub_class.iloc[-1,]["second_category"], atype])
            self.chart_class_column("kegg_annot.",  geneset, source, None, class_source,  None,  "geneset.annot_kegg_stat.bar.json")

    def denovo_chart_geneset_enrich_go(self, go_enrich_table, geneset_list=None, geneset_name="geneset",top=20,p_thre=0.5):
        a = pd.read_table(go_enrich_table, header=0)
        a = a[a["p_corrected"] <= p_thre]
        a['k'] = range(a.shape[0])
        a = a.sort_values(by=["p_corrected", 'k'], ascending=True)[:top]
        # a = a.sort_values(by=["p_uncorrected"], ascending=[True])[:top]
        #a = a[:top]
        # 反转
        # a = a[top-1::-1]
        line_data_value_bar = [int(x.split("/")[0]) for x in a["ratio_in_study"]]
        data_bar = a["neg_log10p_corrected"].tolist()
        category_bar = a["discription"].tolist()
        self.denovo_chart_bar_and_line("go_enrich", ".gene_set",  line_data_value_bar,  data_bar, category_bar, geneset_name, "geneset.enrich_go.bar_line.json")

        source_bar = zip(a['go_type'].tolist(), a['discription'].tolist(), a['study_count'].tolist(), a['neg_log10p_corrected'].tolist())
        source_bar = [list(i) for i in source_bar]
        category = list(set(a['go_type'].tolist()))
        title = geneset_name
        self.denovo_chart_bar("go_enrich", ".gene_set",  source_bar, category, title, "go_enrich_go_bar.json")
        # b = a.sort_values(by=["go_type", "p_uncorrected"], ascending=[False, False])
        # source_bar = [
        #     ["item"] + list(b["discription"]),
        #     ["series"] + list(b["neg_log10p_corrected"]),
        #     ["category"] + list(b["go_type"])
        # ]
        # self.chart_bar("go_enrich", ".gene_set",  source_bar, geneset_name, "geneset.enrich_go.bar.json")

        data_buble1 = list()
        for rec in a.to_dict('records'):
            data_buble1.append([rec["discription"], rec["study_count"], rec["pop_count"], rec["p_corrected"]])

        self.denovo_chart_buble1("go_enrich", ".gene_set",  data_buble1, geneset_name, "go_enrich_bubble1.json")
        data_buble2 = [[], [], []]
        # legend = {'BP': 'Biological Process', 'CC': 'Cellular Component', 'MF': 'Molecular Function'}
        legend_list = list()
        for atype, sub_class in a.groupby("go_type"):
            data_sub = [{'y': rec["neg_log10p_corrected"], 'x': rec["enrich_factor"], 'desc': rec["discription"],
                         'name': rec["go_id"], 'size': rec["study_count"]} for rec in sub_class.to_dict("records")]
            if atype == 'BP':
                data_buble2[0] = data_sub
            if atype == 'CC':
                data_buble2[1] = data_sub
            if atype == 'MF':
                data_buble2[2] = data_sub
            # data_buble2.append(data_sub)
            # legend_list.append({atype: legend[atype]})

        self.denovo_chart_buble2("go_enrich", ".gene_set",  data_buble2, geneset_name, "geneset.enrich_go.buble2.json")

    def denovo_chart_geneset_enrich_kegg(self, kegg_enrich_table, geneset_list=None, geneset_name="geneset", top=20,
                                  p_thre=0.5):
        a = pd.read_table(kegg_enrich_table, header=0)
        a = a[a["Corrected P-Value"] <= p_thre]
        a = a.sort_values(by=["P-Value"], ascending=[True])[:top]
        # 判断最小 log10(p)
        if len([x for x in a["Corrected P-Value"] if x > 0]) > 0:
            pvalues_min = min([x for x in a["Corrected P-Value"] if x > 0]) / 10
        else:
            pvalues_min = 0.0001

        pvalues_min = - math.log10(pvalues_min)
        log10x = [-math.log10(x) if x > 0 else pvalues_min for x in a["Corrected P-Value"]]
        a["neg_log10p_corrected"] = log10x

        a["enrich_factor"] = a["Ratio_in_study"].map(lambda x: float(x.split("/")[0])) / a["Ratio_in_pop"].map(
            lambda x: float(x.split("/")[0]))
        a["study_count"] = a["Ratio_in_study"].map(lambda x: int(x.split("/")[0]))
        a["pop_count"] = a["Ratio_in_pop"].map(lambda x: int(x.split("/")[0]))
        # a = a[top - 1::-1]
        line_data_value_bar = [int(x.split("/")[0]) for x in a["Ratio_in_study"]]
        data_bar = a["neg_log10p_corrected"].tolist()
        category_bar = a["Term"].tolist()
        self.denovo_chart_bar_and_line("kegg_enrich", ".gene_set",  line_data_value_bar,  data_bar, category_bar, geneset_name, "geneset.enrich_go.bar_line.json")

        source_bar = zip(a['typeI'].tolist(), a['Term'].tolist(), a['#Study_num'].tolist(), a['neg_log10p_corrected'].tolist())
        source_bar = [list(i) for i in source_bar]
        category = list(set(a['typeI'].tolist()))
        title = geneset_name
        color = {"Metabolism":"#339933","Environmental Information Processing":"#FF9800","Human Diseases":"#0099FF","Organismal Systems":"#FF00FF","Genetic Information Processing":"#FF0000","Cellular Processes":"#FFFF00","Drug Development":"#666666"}
        self.denovo_chart_bar("kegg_enrich", ".gene_set",  source_bar, category, title, "go_enrich_go_bar.json")

        # kegg_order = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing', 'Cellular Processes', 'Organismal Systems', 'Human Diseases', 'Drug Development']
        # a['order']=[kegg_order.index(x) for x in a['typeI']]
        # b = a.sort_values(by=["order", "neg_log10p_corrected"], ascending=[True, True])
        # source_bar = [
        #     ["item"] + list(b["Term"]),
        #     ["series"] + list(b["neg_log10p_corrected"]),
        #     ["category"] + list(b["typeI"])
        # ]
        # self.chart_bar("kegg_enrich", ".gene_set",  source_bar, geneset_name, "geneset.enrich_go.bar.json")
        # data_buble1 = list()
        # for rec in a.to_dict('records'):
        #     data_buble1.append([rec["Term"], rec["study_count"], rec["enrich_factor"], rec["Corrected P-Value"]])
        #
        # self.denovo_chart_buble1("kegg_enrich", ".gene_set",  data_buble1, geneset_name, "geneset.enrich_go.buble1.json")

        data_buble1 = list()
        for rec in a.to_dict('records'):
            data_buble1.append([rec["Term"], rec["#Study_num"], rec["pop_count"], rec["Corrected P-Value"]])

        self.denovo_chart_buble1("kegg_enrich", ".gene_set",  data_buble1, geneset_name, "go_enrich_bubble1.json")

        data_buble2 = list()
        legend = {'Environmental Information Processing': 'EIP', 'Genetic Information Processing':'GIP', 'Metabolism': 'M',
                  "Cellular Processes": "CP", 'Organismal Systems': "OS", 'Human Diseases': "HD"
                  }
        legend_list = list()
        color_list = list()
        for atype, sub_class in a.groupby("typeI"):
            data_sub = [{'y': rec["neg_log10p_corrected"], 'x': rec["enrich_factor"], 'desc': rec["Term"],
                         'name': rec["ID"], 'size': rec["study_count"]} for rec in sub_class.to_dict("records")]
            data_buble2.append(data_sub)
            legend_list.append({'name': legend[atype], "desc": atype})
            color_list.append(color[atype])

        self.denovo_chart_buble2("kegg_enrich", ".gene_set",  data_buble2, geneset_name, "kegg_enrich_bubble2.json", legend=legend_list, color=color_list)


    def chart_geneset_enrich_go(self, go_enrich_table, geneset_list=None, geneset_name="geneset",top=20,p_thre=0.5):
        a = pd.read_table(go_enrich_table, header=0)
        a = a[a["p_corrected"] <= p_thre]
        # a = a.sort_values(by=["p_uncorrected"], ascending=[True])[:top]
        a = a[:top]
        # 反转
        a = a[top-1::-1]
        source_bar = [
            ["item"] + list(a["discription"]),
            ["series"] + list(a["neg_log10p_corrected"]),
            ["category"] + ["-log10(Padjust)"] * len(a)
        ]

        source_line = [
            {
                "category": "Number",
                "data": [int(x.split("/")[0]) for x in a["ratio_in_study"]],
                "name": "Number",
            }
        ]

        categories_line = list(a["discription"])
        self.chart_bar_and_line("go_enrich", ".gene_set",  source_bar,  source_line, categories_line, geneset_name, "geneset.enrich_go.bar_line.json")

        b = a.sort_values(by=["go_type", "p_uncorrected"], ascending=[False, False])
        source_bar = [
            ["item"] + list(b["discription"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["go_type"])
        ]
        self.chart_bar("go_enrich", ".gene_set",  source_bar, geneset_name, "geneset.enrich_go.bar.json")

        source_buble =  [["x", "y", "size", "fdr"]]
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["discription"],
                rec["study_count"],
                rec["p_corrected"],
            ])

        self.chart_buble("go_enrich", ".gene_set",  source_buble, geneset_name, "geneset.enrich_go.buble.json")

        source_buble_list = list()
        for atype, sub_class in a.groupby("go_type"):
            source =  [["x", "y", "size", "categories", "ID", "Description"]]
            for rec in sub_class.to_dict("records"):
                source.append([
                    rec["enrich_factor"],
                    rec["neg_log10p_corrected"],
                    rec["study_count"],
                    rec["go_type"],
                    rec["go_id"],
                    rec["discription"]
                ])
            source_buble_list.append(source)
        self.chart_buble2("go_enrich", ".gene_set",  source_buble_list, geneset_name, "geneset.enrich_go.buble2.json")

    def chart_geneset_enrich_kegg(self, kegg_enrich_table, geneset_list=None, geneset_name="geneset", top=20,p_thre=0.5):
        a = pd.read_table(kegg_enrich_table, header=0)
        a = a.sort_values(by=["P-Value"], ascending=[True])[:top]
        a = a[a["Corrected P-Value"] <= p_thre]
        # 判断最小 log10(p)
        if len([x for x in a["Corrected P-Value"] if x > 0]) > 0:
            pvalues_min = min([x for x in a["Corrected P-Value"] if x > 0])/10
        else:
            pvalues_min = 0.0001

        pvalues_min = - math.log10(pvalues_min)
        log10x = [-math.log10(x) if x>0 else pvalues_min for x in a["Corrected P-Value"]]
        a["neg_log10p_corrected"] = log10x

        a["enrich_factor"] = a["Ratio_in_study"].map(lambda x: float(x.split("/")[0])) / a["Ratio_in_pop"].map(lambda x: float(x.split("/")[0]))
        a["study_count"] = a["Ratio_in_study"].map(lambda x: int(x.split("/")[0]))

        a = a[top-1::-1]

        source_bar = [
            ["item"] + list(a["Term"]),
            ["series"] + list(a["neg_log10p_corrected"]),
            ["category"] + ["-log10(Padjust)"] * len(a)
        ]

        source_line = [
            {
                "category": "Number",
                "data": list(a["study_count"]),
                "name": "Number",
            }
        ]

        categories_line = list(a["Term"])
        self.chart_bar_and_line("kegg_enrich", ".gene_set",  source_bar,  source_line, categories_line, geneset_name, "geneset.enrich_kegg.bar_line.json")
        kegg_order = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing', 'Cellular Processes', 'Organismal Systems', 'Human Diseases', 'Drug Development']
        a['order']=[kegg_order.index(x) for x in a['typeI']]
        b = a.sort_values(by=["order", "neg_log10p_corrected"], ascending=[True, True])
        source_bar = [
            ["item"] + list(b["Term"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["typeI"])
        ]
        self.chart_bar("kegg_enrich", ".gene_set",  source_bar, geneset_name, "geneset.enrich_kegg.bar.json")

        source_buble =  [["x", "y", "size", "fdr"]]
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["Term"],
                rec["study_count"],
                rec["Corrected P-Value"]
            ])

        self.chart_buble("kegg_enrich", ".gene_set",  source_buble, geneset_name, "geneset.enrich_kegg.buble.json")

        source_buble_list = list()
        for atype, sub_class in a.groupby("typeI"):
            source =  [["x", "y", "size", "categories", "ID", "Description"]]
            for rec in sub_class.to_dict("records"):
                source.append([
                    rec["enrich_factor"],
                    rec["neg_log10p_corrected"],
                    rec["study_count"],
                    rec["typeI"],
                    rec["ID"],
                    rec["Term"]
                ])
            source_buble_list.append(source)
        self.chart_buble2("kegg_enrich", ".gene_set",  source_buble_list, geneset_name, "geneset.enrich_kegg.buble2.json")

    def chart_geneset_enrich_circ(self, circ_table, circ_zscore_input):
        circ_zscore = pd.read_table(circ_zscore_input, header=None)
        term_ids = list(circ_zscore[0])
        term_des = list(circ_zscore[1])
        term_zscores = list(circ_zscore[2])
        
        circ = pd.read_table(circ_table, header=0)
        columns=list(circ.columns)
        columns[0]='seq_id'
        columns[-2]='log2fc'
        columns[-1]='gene_name'
        circ.columns = columns

        source_circ = list()
        for rec in circ.to_dict('records'):
            source_circ.append([rec["gene_name"]] + [rec[term] for term in term_ids] + [rec["log2fc"]])

        source_class = [[t, z] for t,z in zip(term_des, term_zscores)]
        if term_ids[0].startswith("GO"):
            title = "GO term"
        else:
            title = "KEGG term"

        self.chart_circ("geneset_circ", ".enrich", source_circ, source_class, title, "geneset.enrich.circ.json")

    def denovo_chart_geneset_enrich_circ(self, circ_table, circ_zscore_input):
        circ_zscore = pd.read_table(circ_zscore_input, header=None)
        term_ids = list(circ_zscore[0])
        term_des = list(circ_zscore[1])
        term_zscores = list(circ_zscore[2])
        header = [ [term_des[i], term_zscores[i]] for i in range(0, len(term_ids))]

        circ = pd.read_table(circ_table, header=0)
        columns = list(circ.columns)
        columns[0] = 'seq_id'
        columns[-2] = 'log2fc'
        columns[-1] = 'gene_name'
        circ.columns = columns

        source_circ = list()
        for rec in circ.to_dict('records'):
            source_circ.append([rec["gene_name"]] + [rec[term] for term in term_ids] + [rec["log2fc"]])

        if term_ids[0].startswith("GO"):
            title = "GO term"
        else:
            title = "KEGG term"

        self.denovo_chart_circ("geneset_circ", ".enrich", source_circ, header, title, "enrich_circ.json")

    def chart_geneset_enrich_dag(self, dag_table):
        dag_df = pd.read_table(dag_table, header=0)
        dag_df = dag_df.fillna("")

        row_dict_list = dag_df.to_dict('records')
        source = [["child", "parent", "relation", "category", "name", "value", "ratio"]]
        for rec in dag_df.to_dict('records'):
            source.append([
                rec["child"],
                rec["parent"],
                rec["relation"],
                rec["choosed"],
                rec["detail"],
                rec["p_value"],
                None
            ])

        self.chart_dag("geneset_dag", ".dag", source, "geneset.enrich.dag.json")

    def chart_geneset_corr_net(self, net_json):
        with open(net_json, 'r') as f:
            net = json.loads(f.read())
        self.chart_network("geneset_corr", ".net", net, "geneset.corr.net.json")

    def chart_json_batch(self, chart_json):
        if "qc_file" in chart_json:
            qc_files = [chart_json["qc_file"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.chart_raw_qc(chart_json["samples"], qc_files)

        if "align_satu_r" in chart_json and "align_satu_p" in chart_json:
            align_satu_r_files = [chart_json["align_satu_r"].format(sample_name=sample) for sample in chart_json["samples"]]
            align_satu_p_files = [chart_json["align_satu_p"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.chart_satu(chart_json["samples"], align_satu_r_files, align_satu_p_files)

        if "align_coverage" in chart_json:
            align_coverage_files = [chart_json["align_coverage"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.chart_coverage(chart_json["samples"], align_coverage_files)


        if "align_pos" in chart_json:
            align_pos_files = [chart_json["align_pos"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.chart_readpos(chart_json["samples"], align_pos_files)

        if "align_chr" in chart_json:
            align_chr_files = [chart_json["align_chr"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.chart_readchr(chart_json["samples"], align_chr_files)

        if "assemble_step" in chart_json:
            self.chart_assemble(chart_json["assemble_step"])

        if "assemble_new" in chart_json:
            self.chart_assemble_new(chart_json["assemble_new"])

        annot_stat = chart_json["annot_stat"]
        gene_exp = chart_json["gene_exp_all"]
        trans_exp = chart_json["tran_exp_all"]
        venn_dir = os.path.dirname(chart_json["annot_stat"])
        self.chart_annotation_stat("", gene_exp, trans_exp, annot_stat, venn_dir)

        gene_exp = chart_json["gene_exp_ref"]
        tran_exp = chart_json["tran_exp_ref"]
        group_dict = chart_json["group_dict"]
        self.chart_exp_dis("", gene_exp, tran_exp, group_dict)

        self.generate_html_sh()

if __name__ == '__main__':
    a = ChartGeneset()

    # a.chart_geneset_class_cog(cog_class_table=sys.argv[1], geneset_list=None)
    # a.chart_geneset_class_go(go_class_table=sys.argv[1], geneset_list=None)
    # a.chart_splice_diff_stat(splice_diff, splice_psi)
    # a.chart_geneset_class_kegg(kegg_class_table=sys.argv[1])
    # a.chart_geneset_enrich_go(go_enrich_table=sys.argv[1])
    # a.chart_geneset_enrich_kegg(kegg_enrich_table=sys.argv[1])
    # a.chart_geneset_entich_circ(sys.argv[1], sys.argv[2])
    # a.chart_geneset_entich_dag(sys.argv[1])
    # a.chart_geneset_venn("5fb3219b17b2bf1307717d2d,5fb321b017b2bf1307717d2f")

    group_dict = {"A":["A1","A2","A3"], "B": ["B1","B2","B3"], "C": ["C1","C2","C3"], "D": ["D1","D2","D3"]}
    subcluster_list = glob.glob("subcluster_*")
    a.chart_geneset_cluster(cluster_exp=sys.argv[1], cluster_tree=sys.argv[2], sample_tree=sys.argv[3], subcluster_list=subcluster_list, group_dict=group_dict, samples_order=["A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3"])
