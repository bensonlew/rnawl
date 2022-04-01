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
from collections import OrderedDict,defaultdict
from chart import Chart
from bson.objectid import ObjectId
import glob


class ChartGeneset(Chart):
    def __init__(self):
        super(ChartGeneset, self).__init__()

    def chart_geneset_venn(self, geneset_ids):
        print(geneset_ids)
        project_type = 'medical_transcriptome'
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

        self.chart_venn("geneset", "", source, "geneset.venn.venn.json")
        self.chart_upset("geneset", "", source, "geneset.venn.upset.json")

    def chart_geneset_cluster(self, cluster_exp, cluster_tree, sample_tree, subcluster_list, group_dict=None, samples_order=None,seq_id2name = None):
        cluster_pd = pd.read_table(cluster_exp, index_col=0, header=0)
        gene_num = cluster_pd.shape[0]
        chart_height = 600 if gene_num <= 600 else 1500

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
                [seq_id2name[line_dict[0]]] + [line_dict[1][s] for s in samples]
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
        gene2group_dict = OrderedDict()
        sub_2gene = OrderedDict()
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
            sub_2gene[c_num] = sub_genes

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

        for sub in sorted(sub_2gene):
            for gene in sub_2gene[sub]:
                gene2group_source.append([gene, gene2group_dict[gene]])
        # for gene in genes:
        #     if gene in gene2group_dict:
        #         gene2group_source.append([gene, gene2group_dict[gene]])
        # # sub_order =sorted(list(set([k[1] for k in gene2group_source])))
        # gene2group_source = sorted(gene2group_source, key = lambda k: k[1])

        self.chart_heat_tree("geneset", ".cluster", corr_heat, sample_tree, gene_tree,  sample2group_source, gene2group_source, "../medical_transcriptome/geneset.cluster.heatmap.json",height=chart_height)

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
            for g, ss in group_dict.items():
                for s in ss:
                    sample2group_source.append([s, g])
        else:
            for s in samples:
                sample2group_source.append([s, s])

        gene2group_source = [["name", "group"]]
        for gene in genes:
            gene2group_source.append([gene, gene])

        self.chart_heat_tree("geneset_gsva", ".cluster", corr_heat, sample_tree, gene_tree, sample2group_source,
                             gene2group_source, "geneset.gsva.cluster.heatmap.json")

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
        up_dict = {'data': [], 'name': 'up'}
        nosig_dict = {'data': [], 'name': 'nosig'}
        down_dict = {'data': [], 'name': 'down'}
        for i in geneset_list:
            if limma_20.loc[i]['padjust'] >= 0.05:
                nosig_dict['data'].append(limma_20.loc[i]['log2fc'])
            else:
                if limma_20.loc[i]['log2fc'] >= 0.2:
                    up_dict['data'].append(limma_20.loc[i]['log2fc'])
                if -0.2 < limma_20.loc[i]['log2fc'] < 0.2:
                    nosig_dict['data'].append(limma_20.loc[i]['log2fc'])
                if limma_20.loc[i]['log2fc'] <= -0.2:
                    down_dict['data'].append(limma_20.loc[i]['log2fc'])
        self.chart_gsva_diff('geneset_gsva_diff.', name, up_dict, nosig_dict, down_dict, geneset_list,
                             'geneset.gsva.diff.json')

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
            geneset_list = [c.split("_COG")[0] for c in b.columns if c.endswith("_COG")]

        source = list()
        for geneset in geneset_list:
            source.append({
                "data": list(b[geneset + "_COG"]),
                "name": geneset
            })
        self.chart_class_column("cog_annot", ".gene_set", source, categories, None, None,
                                "geneset.annot_cog_stat.bar.json")

    def chart_geneset_class_go(self, go_class_table, geneset_list=None, top=20):
        a = pd.read_table(go_class_table, header=0, index_col=0)
        # 原表格多一列处理
        col_names = a.columns
        b = a.iloc[:, :-1]
        b.columns = col_names[1:]
        b["sum"] = sum([b[x] for x in col_names if x.endswith("num")])
        # c = b.sort_values(by=["sum"], ascending=[False])[:top]
        # d = c.sort_index(ascending=False)
        # d["type"] = d.index
        # d = d.sort_values(by=["type", "sum"], ascending=[False, False])
        # categories = list(d["Term"])

        b["k"] = range(b.shape[0])
        c = b.sort_values(by=["sum", "k"], ascending=[False,True])[:top]
        d = c.sort_index(ascending=False)
        d["type"] = d.index
        d = d.sort_values(by=["type", "sum", "k"], ascending=[False, False, True])
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
        for atype, sub_class in d.groupby("type", sort=False):
            class_source.append([sub_class.iloc[0,]["Term"], sub_class.iloc[-1,]["Term"], atype])
        if len(geneset_list) >= 2:
            bac_json = "geneset.annot_go_stat.bar.json"
        else:
            bac_json = "geneset.annot_go_stat2.bar.json"

        self.chart_class_column("go_annot", ".gene_set", source, categories, class_source, None, bac_json)

    def chart_geneset_class_do(self, do_class_table, geneset_list=None, top=20):
        a = pd.read_table(do_class_table, header=0)
        # 原表格多一列处理
        col_names = a.columns
        # b = a.iloc[:, :-1]
        # b.columns = col_names[1:]
        b = a
        b["sum"] = sum([b[x] for x in col_names if x.endswith("numbers")])
        c = b.sort_values(by=["sum"], ascending=[False])[:top]
        d = c.sort_index(ascending=False)
        d["type"] = d.index
        d = d.sort_values(by=["sum"], ascending=[False])
        type_list = list(d['Term Type'])
        # 分类按第一个属于的num排序
        d['k'] = d["Term Type"].map(lambda x: type_list.index(x))
        # d = d.sort_values(by=["k"], ascending=[True])
        d = d.sort_values(by=["k", "sum", "type"], ascending=[True, False, True])
        categories = list(d["DO Name"])
        if geneset_list:
            pass
        else:
            geneset_list = [c.split(" ")[0] for c in a.columns if c.endswith("numbers")]

        source = list()

        print
        geneset_list
        if len(geneset_list) == 1:
            # 单基因集
            geneset = geneset_list[0]
            source = [
                ["item"] + list(d["DO Name"]),
                ["series"] + list(d[geneset + " numbers"]),
                ["category"] + list(d["Term Type"])
            ]
        else:
            for geneset in geneset_list:
                try:
                    percent = [float(p.split("(")[0]) for p in d[geneset + " percent"]]
                except:
                    percent = [float(0) for p in d[geneset + " percent"]]
                source.append({
                    "percent": percent,
                    "data": list(d[geneset + " numbers"]),
                    "name": geneset
                })

        class_source = [["p1", "p2", "categoy"]]
        for atype, sub_class in d.groupby("Term Type", sort=False):
            class_source.append([sub_class.iloc[0,]["DO Name"], sub_class.iloc[-1,]["DO Name"], atype])
        if len(geneset_list) >= 2:
            bac_json = "../medical_transcriptome/geneset.annot_do_stat.bar.json"
        else:
            bac_json = "../medical_transcriptome/geneset.annot_do_stat2.bar.json"

        self.chart_class_column("do_annot", ".gene_set", source, categories, class_source, None, bac_json, reset_margin=False)

    def chart_geneset_class_reactome(self, reactome_class_table, geneset_list=None, top=20):
        a = pd.read_table(reactome_class_table, header=0)
        # 原表格多一列处理
        col_names = a.columns
        b = a.iloc[:, :-1]
        b.columns = col_names[1:]
        b = a
        b["sum"] = sum([b[x] for x in col_names if x.endswith("numbers")])
        c = a.sort_values(by=["sum"], ascending=[False])[:top]
        # d = c.sort_index(ascending=False)
        # d["type"] = d.index
        # d = d.sort_values(by=["sum"], ascending=[False])
        # type_list = list(d['Term Type'])
        # # 分类按第一个属于的num排序
        # d['k'] = d["Term Type"].map(lambda x:type_list.index(x))
        # d = d.
        # sort_values(by=["k"], ascending=[True])
        d = a[a["Category Function description"].isin(list(c["Category Function description"]))]
        categories = list(d["Category Function description"])
        if geneset_list:
            pass
        else:
            geneset_list = [c.split("_numbers")[0] for c in a.columns if c.endswith("numbers")]

        source = list()

        print
        geneset_list
        if len(geneset_list) == 1:
            # 单基因集
            geneset = geneset_list[0]
            source = [
                ["item"] + list(d["Category Function description"]),
                ["series"] + list(d[geneset + "_numbers"]),
                ["category"] + list(d["Category Function description"])
            ]
        else:
            for geneset in geneset_list:
                source.append({
                    "pathway": list(d['Category Function description']),
                    "data": list(d[geneset + "_numbers"]),
                    "name": geneset
                })

        # class_source = [["p1", "p2", "categoy"]]
        # for atype,sub_class in d.groupby("Term Type", sort=False):
        #     class_source.append([sub_class.iloc[0,]["DO Name"], sub_class.iloc[-1,]["DO Name"], atype])
        class_source = None
        if len(geneset_list) >= 2:
            bac_json = "../medical_transcriptome/geneset.annot_reactome_stat.bar.json"
        else:
            bac_json = "../medical_transcriptome/geneset.annot_reactome_stat2.bar.json"

        self.chart_class_column("reactome_annot", ".gene_set", source, categories, class_source, None, bac_json, reset_margin=False)

    def convert_kegg_class_file(self,work_dir,geneset,kegg_class_file,kegg_level_path):
        if not os.path.exists(os.path.join(work_dir,"temporary",geneset)):
            os.makedirs(os.path.join(work_dir,"temporary",geneset))
        work_dir = os.path.join(work_dir,"temporary",geneset)
        stat = pd.read_table(kegg_class_file, header=0)

        # bb = stat
        # col_names = bb.columns
        # bb["sum"] = sum([bb[x] for x in col_names if x.endswith("num")])
        # bb["k"] = range(bb.shape[0])
        # stat =  bb.sort_values(by=["sum", "k"], ascending=[False,True])[:20]
        # stat = stat.drop(["sum","k"],axis=1)
        level = pd.read_table(kegg_level_path, header=0)
        stat_level = pd.merge(stat, level, on='Pathway_id')
        stat_level.to_csv(work_dir + "/" + "stat_level", sep='\t', index=False)

        # 按照kegg官网进行一级分类的排序
        list_custom = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing',
                       'Cellular Processes', 'Organismal Systems', 'Human Diseases',
                       'Drug Development']
        appended_data = []
        first_category = {i for i in stat_level.first_category}
        for i in list_custom:
            if i in first_category:
                data = stat_level.loc[stat_level['first_category'] == i]
                appended_data.append(data)

        appended_data = pd.concat(appended_data)
        appended_data.drop(['graph_id', 'hyperlink', 'graph_png_id'], axis=1, inplace=True)
        appended_data.to_csv(work_dir + "/" + "kegg_annotation_analysis", sep='\t', index=False)
        ls = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep='\t')
        if len(ls.columns) <12 :
            # if not regulate :
            with open(work_dir + "/" + "kegg_annotation_analysis") as f_kaa, open(
                    work_dir + "/" + "kegg_annotation_analysis_new", "w") as fw_kaa:
                header_kaa = f_kaa.readline()
                hkl = header_kaa.strip().split("\t")
                geneko_name1 = hkl[3][:-1] + "ko"

                fw_kaa.write(
                    hkl[0] + "\t" + hkl[1] + "\t" + hkl[2] + "\t" + geneko_name1 + "\t" + hkl[3] + "\t" + hkl[4] + "\t"
                    + hkl[5] + "\t" +
                    hkl[6] + "\t" + hkl[7] + "\t" + hkl[8] + "\t" + hkl[9] + "\t" + hkl[10] + "\n")

                for line in f_kaa:
                    genes_1 = []
                    genes_2 = []
                    line_list = line.strip().split("\t")
                    line_list_3 = line_list[3]
                    name1_genes = line_list_3.split(');')
                    genes_1 += [x.split('(')[0] for x in name1_genes]

                    fw_kaa.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + line_list[3] + "\t" +
                                 ",""".join(genes_1) + "\t" + line_list[4] + "\t" + line_list[5] + "\t" + line_list[
                                     6] + "\t" +
                                 line_list[7] + "\t" + line_list[8] + "\t" + line_list[9] + "\t" + line_list[10] + "\n")
            new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep='\t', header=0)
            new_data.groupby("second_category")
            group_obj = new_data.groupby("second_category")
            groups = group_obj.groups.keys()
            # 做一个基因集
            genesets = new_data.columns[3].split()
            result = defaultdict(dict)
            for each in groups:
                first = new_data.loc[new_data["second_category"] == each]['first_category']
                first = first.to_dict().values()[0]
                for geneset in genesets:
                    group_detail = group_obj.get_group(each)
                    genes = list()
                    for g in group_detail[geneset]:

                        if not pd.isnull(g):
                            tmp = g.split(');')
                            genes += [x.split('(')[0] for x in tmp]
                            # genes = [i for i in genes if genes.count(i) == 1]
                            genes = list(set(genes))
                        else:
                            genes = []
                    result[geneset][each] = [len(genes), first]

            try:
                a = pd.DataFrame(result)
                a.reset_index(inplace=True)
                a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
                a.to_csv(work_dir + "/" + "k", sep='\t', index=False)
                with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                    header = f1.readline()
                    geneset_name1 = header.strip().split("\t")[1]
                    geneset_name1_num = geneset_name1 + "_num"
                    fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + geneset_name1_num + "\n")
                    for line in f1:
                        line_split = line.strip().split("\t")
                        sec = line_split[0]
                        num1 = line_split[1].strip("[]").split(",")[0]
                        first_cate = line_split[1].strip("[]").split(",")[1].strip().strip("'")
                        fw.write(first_cate + "\t" + sec + "\t" + num1 + "\n")
                df_a = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")

                list_custom = ['Metabolism', 'Genetic Information Processing',
                               'Environmental Information Processing',
                               'Cellular Processes', 'Organismal Systems',
                               'Human Diseases',
                               'Drug Development']
                appended_data_new_1 = []
                for i in list_custom:
                    if i in list(df_a.first_category):
                        data = df_a.loc[df_a['first_category'] == i]
                        appended_data_new_1.append(data)

                appended_data_new_1 = pd.concat(appended_data_new_1)
                appended_data_new_1.to_csv(work_dir + "/" + "kegg_statistic", sep='\t', index=False)
            except:
                pass


        else:
            with open(work_dir + "/" + "kegg_annotation_analysis") as f_kaa, open(
                    work_dir + "/" + "kegg_annotation_analysis_new", "w") as fw_kaa:
                header_kaa = f_kaa.readline()
                hkl = header_kaa.strip().split("\t")
                geneko_name1 = hkl[3][:-1] + "ko"
                geneko_name2 = hkl[5][:-1] + "ko"
                fw_kaa.write(hkl[0] + "\t" + hkl[1] + "\t" + hkl[2] + "\t" + geneko_name1 + "\t" + hkl[3] + "\t" + hkl[
                    4] + "\t" +
                             geneko_name2 + "\t" + hkl[5] + "\t" + hkl[6] + "\t" + hkl[7] + "\t" + hkl[8] + "\t" + hkl[
                                 9] + "\t" + hkl[10] +
                             "\t" + hkl[11] + "\t" + hkl[12] + "\n")

                for line in f_kaa:
                    genes_1 = []
                    genes_2 = []
                    line_list = line.strip().split("\t")
                    line_list_3 = line_list[3]
                    name1_genes = line_list_3.split(');')
                    genes_1 += [x.split('(')[0] for x in name1_genes]

                    line_list_5 = line_list[5]
                    name2_genes = line_list_5.split(');')
                    genes_2 += [x.split('(')[0] for x in name2_genes]
                    fw_kaa.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + line_list[3] + "\t" +
                                 ",""".join(genes_1) + "\t" + line_list[4] + "\t" + line_list[5] + "\t" + ",".join(
                        genes_2) + "\t" +
                                 line_list[6] + "\t" + line_list[7] + "\t" + line_list[8] + "\t" + line_list[9] + "\t" +
                                 line_list[10] + "\t" + line_list[11] + "\t" + line_list[12] + "\n")
            new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep='\t', header=0)
            new_data.groupby("second_category")
            group_obj = new_data.groupby("second_category")
            groups = group_obj.groups.keys()
            # 做2个基因集
            genesets = new_data.columns[3], new_data.columns[5]
            result = defaultdict(dict)
            for each in groups:
                first = new_data.loc[new_data["second_category"] == each]['first_category']
                first = first.to_dict().values()[0]
                for geneset in genesets:
                    group_detail = group_obj.get_group(each)
                    genes = list()
                    for g in group_detail[geneset]:

                        # isnull支持的数据类型更多，相比isnan
                        if not pd.isnull(g):
                            tmp = g.split(');')
                            genes += [x.split('(')[0] for x in tmp]
                            # 用set会弹出不知名的错误
                            # genes = [i for i in genes if genes.count(i) == 1]
                        else:
                            pass
                            #  genes = []
                    genes = list(set(genes))
                    result[geneset][each] = [len(genes), first]
            # try:
            a = pd.DataFrame(result)
            a.reset_index(inplace=True)
            a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
            a.to_csv(work_dir + "/" + "k", sep='\t', index=False)
            with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                header = f1.readline()
                geneset_name1 = header.strip().split("\t")[1]
                geneset_name1_num = geneset_name1 + "_num"
                geneset_name2 = header.strip().split("\t")[2]
                geneset_name2_num = geneset_name2 + "_num"
                fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + geneset_name1_num + "\t" + \
                         geneset_name2_num + "\n")
                for line in f1:
                    line_split = line.strip().split("\t")
                    sec = line_split[0]
                    num1 = line_split[1].strip("[]").split(",")[0]
                    num2 = line_split[2].strip("[]").split(",")[0]
                    first_cate = line_split[1].strip("[]").split(",")[1].strip().strip("'")
                    fw.write(first_cate + "\t" + sec + "\t" + num1 + "\t" + num2 + "\n")
            df_a = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")

            list_custom = ['Metabolism', 'Genetic Information Processing',
                           'Environmental Information Processing',
                           'Cellular Processes', 'Organismal Systems',
                           'Human Diseases',
                           'Drug Development']
            appended_data_new_2 = []
            for i in list_custom:
                if i in list(df_a.first_category):
                    data = df_a.loc[df_a['first_category'] == i]
                    appended_data_new_2.append(data)
            appended_data_new_2 = pd.concat(appended_data_new_2)
            appended_data_new_2.to_csv(work_dir + "/" +
                                       "kegg_statistic", sep='\t', index=False)
        return     work_dir + "/" + "kegg_statistic"

    def chart_geneset_class_kegg(self, kegg_class_table, geneset_list=None, top=20):
        a = pd.read_table(kegg_class_table, header=0)
        # a = a[:top]

        geneset_list = [c.split("_num")[0] for c in a.columns if c.endswith("num")]
        for geneset in geneset_list:
            a_choose = a[a[geneset + "_num"] > 0]
            a_choose["id"] = a_choose.index
            i_choose = a_choose.sort_values(by=[geneset + "_num", "id"], ascending=[False, True])[:top]
            # kegg_order = list(i_choose.drop_duplicates("first_category")["first_category"])
            kegg_order = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing',
                          'Cellular Processes', 'Organismal Systems', 'Human Diseases', 'Drug Development']
            # categories = list(i_choose["first_category"])
            i_choose['order'] = [kegg_order.index(x) for x in i_choose['first_category']]
            f_choose = i_choose.sort_values(by=["order", geneset + "_num"], ascending=[True, False])
            # a_choose = a[a[geneset + "_num"] > 0]
            source = [
                ["item"] + list(f_choose["second_category"]),
                ["series"] + list(f_choose[geneset + "_num"]),
                ["category"] + list(f_choose["first_category"])
            ]

            class_source = ["p1", "p2", "category"]
            for atype, sub_class in a.groupby("first_category"):
                class_source.append(
                    [sub_class.iloc[0,]["second_category"], sub_class.iloc[-1,]["second_category"], atype])
            self.chart_class_column("kegg_annot.", geneset, source, None, class_source, None,
                                    "geneset.annot_kegg_stat.bar.json")

    def chart_geneset_enrich_go(self, go_enrich_table, geneset_list=None, geneset_name="geneset", top=20, p_thre=0.5):
        a = pd.read_table(go_enrich_table, header=0)
        a = a[a["p_corrected"] <= p_thre]
        a["k"] = range(a.shape[0])
        # a = a.sort_values(by=["p_uncorrected"], ascending=[True])[:top]
        a = a.sort_values(by=["p_corrected", "k"], ascending=[True, True])[:top]
        # 反转
        a = a[top - 1::-1]
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
        self.chart_bar_and_line("go_enrich", ".gene_set", source_bar, source_line, categories_line, geneset_name,
                                "geneset.enrich_go.bar_line.json")

        b = a.sort_values(by=["go_type", "p_uncorrected"], ascending=[False, False])
        source_bar = [
            ["item"] + list(b["discription"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["go_type"])
        ]
        self.chart_bar("go_enrich", ".gene_set", source_bar, geneset_name, "geneset.enrich_go.bar.json")

        source_buble = [["x", "y", "size", "fdr"]]
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["discription"],
                rec["study_count"],
                rec["p_corrected"],
            ])

        self.chart_buble("go_enrich", ".gene_set", source_buble, geneset_name, "geneset.enrich_go.buble.json")

        source_buble_list = list()
        for atype, sub_class in a.groupby("go_type"):
            source = [["x", "y", "size", "categories", "ID", "Description"]]
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
        self.chart_buble2("go_enrich", ".gene_set", source_buble_list, geneset_name, "geneset.enrich_go.buble2.json")

    def chart_geneset_enrich_kegg(self, kegg_enrich_table, geneset_list=None, geneset_name="geneset", top=20,
                                  p_thre=0.5):
        a = pd.read_table(kegg_enrich_table, header=0)
        # a = a.sort_values(by=["P-Value"], ascending=[True])[:top]
        a["k"] = range(a.shape[0])
        a = a.sort_values(by=["Corrected P-Value","k"], ascending=[True,True])[:top]
        a = a[a["Corrected P-Value"] <= p_thre]
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

        a = a[top - 1::-1]

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
        self.chart_bar_and_line("kegg_enrich", ".gene_set", source_bar, source_line, categories_line, geneset_name,
                                "geneset.enrich_kegg.bar_line.json")


        # kegg_order = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing',
        #               'Cellular Processes', 'Organismal Systems', 'Human Diseases', 'Drug Development']
        # a['order'] = [kegg_order.index(x) for x in a['typeI']]
        # b = a.sort_values(by=["order", "neg_log10p_corrected"], ascending=[True, True])
        # source_bar = [
        #     ["item"] + list(b["Term"]),
        #     ["series"] + list(b["neg_log10p_corrected"]),
        #     ["category"] + list(b["typeI"])
        # ]
        # 不用管这段代码干嘛的,有问题去看有参原版,这版为了适应医学的产品端zz需求弄的
        t = a[top - 1::-1]
        kegg_order = list(t.drop_duplicates("typeI")["typeI"])

        def Reverse(lst):
            return [ele for ele in reversed(lst)]

        kegg_order = Reverse(kegg_order)
        t['order'] = [kegg_order.index(x) for x in t['typeI']]

        # kegg_order =list(a.drop_duplicates("typeI")["typeI"])
        # a['order']=[kegg_order.index(x) for x in a['typeI']]
        b = t.sort_values(by=["order", "neg_log10p_corrected"], ascending=[True, True])
        kegg_relate = {
            'Metabolism': "M",
            'Genetic Information Processing': "GIP",
            'Environmental Information Processing': "EIP",
            'Cellular Processes': "CP",
            'Organismal Systems': "OS",
            'Human Diseases': "HD",
            "Drug Development": "DD"
        }
        # b = a.sort_values(by=["P-Value"], ascending=[True])[:top]
        b["typeI"] = b.apply(lambda x: kegg_relate[x["typeI"]], axis=1)

        source_bar = [
            ["item"] + list(b["Term"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["typeI"])
        ]
        self.chart_bar("kegg_enrich", ".gene_set", source_bar, geneset_name, "geneset.enrich_kegg.bar.json")

        source_buble = [["x", "y", "size", "fdr"]]
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["Term"],
                rec["study_count"],
                rec["Corrected P-Value"]
            ])

        self.chart_buble("kegg_enrich", ".gene_set", source_buble, geneset_name, "geneset.enrich_kegg.buble.json")

        source_buble_list = list()
        for atype, sub_class in a.groupby("typeI"):
            source = [["x", "y", "size", "categories", "ID", "Description"]]
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
        self.chart_buble2("kegg_enrich", ".gene_set", source_buble_list, geneset_name,
                          "geneset.enrich_kegg.buble2.json")

    def chart_geneset_enrich_do(self, do_enrich_table, geneset_list=None, geneset_name="geneset", top=20, p_thre=0.5):
        a = pd.read_table(do_enrich_table, header=0)
        a = a[a["Padjust"] <= p_thre]
        a = a.sort_values(by=["Padjust"], ascending=[True])[:top]

        # 判断最小 log10(p)
        if len([x for x in a["Padjust"] if x > 0]) > 0:
            pvalues_min = min([x for x in a["Padjust"] if x > 0]) / 10
        else:
            pvalues_min = 0.0001

        pvalues_min = - math.log10(pvalues_min)
        log10x = [-math.log10(x) if x > 0 else pvalues_min for x in a["Padjust"]]
        a["neg_log10p_corrected"] = log10x
        a["enrich_factor"] = a.apply(
            lambda x: float(x["Ratio_in_study"].split("/")[0]) / float(x["Ratio_in_pop"].split("/")[0]), axis=1)
        a["study_count"] = a.apply(lambda x: int(x["Ratio_in_study"].split("/")[0]), axis=1)

        # 反转
        a = a[top - 1::-1]
        source_bar = [
            ["item"] + list(a["DO name"]),
            ["series"] + list(a["neg_log10p_corrected"]),
            ["category"] + ["-log10(Padjust)"] * len(a)
        ]

        source_line = [
            {
                "category": "Number",
                "data": [int(x.split("/")[0]) for x in a["Ratio_in_study"]],
                "name": "Number",
            }
        ]

        categories_line = list(a["DO name"])
        self.chart_bar_and_line("do_enrich", ".gene_set", source_bar, source_line, categories_line, geneset_name,
                                "../medical_transcriptome/geneset.enrich_do.bar_line.json")

        b = a.sort_values(by=["Term_type", "Padjust"], ascending=[False, False])
        source_bar = [
            ["item"] + list(b["DO name"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["Term_type"])
        ]
        self.chart_bar("do_enrich", ".gene_set", source_bar, geneset_name,
                       "../medical_transcriptome/geneset.enrich_do.bar.json")

        source_buble = [["x", "y", "size", "fdr"]]
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["DO name"],
                rec["study_count"],
                rec["Padjust"],
            ])

        self.chart_buble("do_enrich", ".gene_set", source_buble, geneset_name,
                         "../medical_transcriptome/geneset.enrich_do.buble.json")

        source_buble_list = list()
        for atype, sub_class in a.groupby("Term_type"):
            source = [["x", "y", "size", "categories", "ID", "Description"]]
            for rec in sub_class.to_dict("records"):
                source.append([
                    rec["enrich_factor"],
                    rec["neg_log10p_corrected"],
                    rec["study_count"],
                    rec["Term_type"],
                    rec["DO_ID"],
                    rec["DO name"]
                ])
            source_buble_list.append(source)
        self.chart_buble2("do_enrich", ".gene_set", source_buble_list, geneset_name,
                          "../medical_transcriptome/geneset.enrich_do.buble2.json")

    def chart_geneset_enrich_reactome(self, reactome_enrich_table, geneset_list=None, geneset_name="geneset", top=20,
                                      p_thre=0.5):
        a = pd.read_table(reactome_enrich_table, header=0)
        a = a[a["Padjust"] <= p_thre]
        a = a.sort_values(by=["Padjust"], ascending=[True])[:top]

        # 判断最小 log10(p)
        if len([x for x in a["Padjust"] if x > 0]) > 0:
            pvalues_min = min([x for x in a["Padjust"] if x > 0]) / 10
        else:
            pvalues_min = 0.0001

        pvalues_min = - math.log10(pvalues_min)
        log10x = [-math.log10(x) if x > 0 else pvalues_min for x in a["Padjust"]]
        a["neg_log10p_corrected"] = log10x
        if len(a) == 0:
            a["enrich_factor"] = []
            a["study_count"] = []
        else:
            a["enrich_factor"] = a.apply(
                lambda x: float(x["Ratio_in_study"].split("/")[0]) / float(x["Ratio_in_pop"].split("/")[0]), axis=1)
            a["study_count"] = a.apply(lambda x: int(x["Ratio_in_study"].split("/")[0]), axis=1)

        # 反转
        a = a[top - 1::-1]
        source_bar = [
            ["item"] + list(a["Description"]),
            ["series"] + list(a["neg_log10p_corrected"]),
            ["category"] + ["-log10(Padjust)"] * len(a)
        ]

        source_line = [
            {
                "category": "Number",
                "data": [int(x.split("/")[0]) for x in a["Ratio_in_study"]],
                "name": "Number",
            }
        ]

        categories_line = list(a["Description"])
        self.chart_bar_and_line("reactome_enrich", ".gene_set", source_bar, source_line, categories_line, geneset_name,
                                "../medical_transcriptome/geneset.enrich_reactome.bar_line.json")

        b = a.sort_values(by=["category", "Padjust"], ascending=[False, False])
        source_bar = [
            ["item"] + list(b["Description"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["category"])
        ]
        self.chart_bar("reactome_enrich", ".gene_set", source_bar, geneset_name,
                       "../medical_transcriptome/geneset.enrich_reactome.bar.json")

        source_buble = [["x", "y", "size", "fdr"]]
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["Description"],
                rec["study_count"],
                rec["Padjust"],
            ])

        self.chart_buble("reactome_enrich", ".gene_set", source_buble, geneset_name,
                         "../medical_transcriptome/geneset.enrich_reactome.buble.json")

        source_buble_list = list()
        for atype, sub_class in a.groupby("category"):
            source = [["x", "y", "size", "categories", "ID", "Description"]]
            for rec in sub_class.to_dict("records"):
                source.append([
                    rec["enrich_factor"],
                    rec["neg_log10p_corrected"],
                    rec["study_count"],
                    rec["category"],
                    rec["Pathway ID"],
                    rec["Description"]
                ])
            source_buble_list.append(source)
        self.chart_buble2("reactome_enrich", ".gene_set", source_buble_list, geneset_name,
                          "../medical_transcriptome/geneset.enrich_reactome.buble2.json")

    def chart_geneset_enrich_disgenet(self, disgenet_enrich_table, geneset_name="geneset", top=20):
        a = pd.read_table(disgenet_enrich_table, header=0)

        # 判断最小 log10(p)
        if len([x for x in a["p.adjust"] if x > 0]) > 0:
            pvalues_min = min([x for x in a["p.adjust"] if x > 0]) / 10
        else:
            pvalues_min = 0.0001

        pvalues_min = - math.log10(pvalues_min)
        log10x = [-math.log10(x) if x > 0 else pvalues_min for x in a["p.adjust"]]
        a["neg_log10p_corrected"] = log10x

        a = a[top - 1::-1]
        source_bar = [
            ["item"] + list(a["Description"]),
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

        categories_line = list(a["Description"])
        self.chart_bar_and_line("disgenet_enrich", ".gene_set", source_bar, source_line, categories_line, geneset_name,
                                "../medical_transcriptome/geneset.enrich_disgenet.bar_line.json")

        b = a.sort_values(by=["Description", "p.adjust"], ascending=[False, False])
        source_bar = [
            ["item"] + list(b["Description"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["Description"])
        ]
        self.chart_bar("disgenet_enrich", ".gene_set", source_bar, geneset_name,
                       "../medical_transcriptome/geneset.enrich_disgenet.bar.json")

        source_buble = [["x", "y", "size", "fdr"]]
        a["enrich_factor"] = a.apply(
            lambda x: float(x["ratio_in_study"].split("/")[0]) / float(x["ratio_in_pop"].split("/")[0]), axis=1)
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["Description"],
                rec["Count"],
                rec["p.adjust"],
            ])

        self.chart_buble("disgenet_enrich", ".gene_set", source_buble, geneset_name,
                         "../medical_transcriptome/geneset.enrich_disgenet.bubble.json")

        source_buble_list = list()
        a['category'] = 'disease'
        for atype, sub_class in a.groupby("category"):
            source = [["x", "y", "size", "categories", "ID", "Description"]]
            for rec in sub_class.to_dict("records"):
                source.append([
                    rec["enrich_factor"],
                    rec["neg_log10p_corrected"],
                    rec["Count"],
                    rec['category'],
                    rec["ID"],
                    rec["Description"]
                ])
            source_buble_list.append(source)
        self.chart_buble2("disgenet_enrich", ".gene_set", source_buble_list, geneset_name,
                          "../medical_transcriptome/geneset.enrich_disgenet.bubble2.json")

    def chart_geneset_enrich_circ(self, circ_table, circ_zscore_input):
        circ_zscore = pd.read_table(circ_zscore_input, header=None)
        term_ids = list(circ_zscore[0])
        term_des = list(circ_zscore[1])
        term_zscores = list(circ_zscore[2])

        circ = pd.read_table(circ_table, header=0)
        columns = list(circ.columns)
        columns[0] = 'seq_id'
        columns[-2] = 'log2fc'
        columns[-1] = 'gene_name'
        circ.columns = columns

        source_circ = list()
        for rec in circ.to_dict('records'):
            source_circ.append([rec["gene_name"]] + [rec[term] for term in term_ids] + [rec["log2fc"]])

        source_class = [[t, z] for t, z in zip(term_des, term_zscores)]
        if term_ids[0].startswith("GO"):
            title = "GO term"
        elif term_ids[0].startswith("DOID"):
            title = "DO term"
        elif term_ids[0].startswith("R-"):
            title = "Reactome term"
        elif term_ids[0].startswith("C"):
            title = "DisGeNET term"
        else:
            title = "KEGG term"

        self.chart_circ("geneset_circ", ".enrich", source_circ, source_class, title, "geneset.enrich.circ.json")

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
            align_satu_r_files = [chart_json["align_satu_r"].format(sample_name=sample) for sample in
                                  chart_json["samples"]]
            align_satu_p_files = [chart_json["align_satu_p"].format(sample_name=sample) for sample in
                                  chart_json["samples"]]
            self.chart_satu(chart_json["samples"], align_satu_r_files, align_satu_p_files)

        if "align_coverage" in chart_json:
            align_coverage_files = [chart_json["align_coverage"].format(sample_name=sample) for sample in
                                    chart_json["samples"]]
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

    # group_dict = {"A":["A1","A2","A3"], "B": ["B1","B2","B3"], "C": ["C1","C2","C3"], "D": ["D1","D2","D3"]}
    # subcluster_list = glob.glob("subcluster_*")
    # a.chart_geneset_cluster(cluster_exp=sys.argv[1], cluster_tree=sys.argv[2], sample_tree=sys.argv[3], subcluster_list=subcluster_list, group_dict=group_dict, samples_order=["A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3"])
    # a.chart_geneset_class_reactome(reactome_class_table=sys.argv[1], geneset_list=None)
    # a.chart_geneset_enrich_do(do_enrich_table=sys.argv[1])
    a.chart_geneset_enrich_reactome(reactome_enrich_table=sys.argv[1])

