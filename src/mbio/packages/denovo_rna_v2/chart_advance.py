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
from mbio.packages.rna.dendrogram2newick import convert2
from matplotlib import colors

class ChartAdvance_denovo(Chart):
    def __init__(self):
        super(ChartAdvance_denovo, self).__init__()

    def chart_wgcna_sample_tree(self, sample_tree_table, group_dict=None):
        linkage_pd = pd.read_table(sample_tree_table, header=None)
        linkage_list = [list(x) for x in linkage_pd.as_matrix()]
        ordered_leaf = pd.read_table(sample_tree_table[:-3] + 'order.txt', header=None)[0]
        ordered_seqs = [x.strip() for x in open(sample_tree_table[:-3] + 'labels.txt')]
        newick_tree = convert2(sample_tree_table, ordered_leaf, ordered_seqs)
        source = list()

        sample2group_source = [["name", "group"]]
        if group_dict:
            for g,ss in group_dict.items():
                for s in ss:
                    sample2group_source.append([s, g])
        else:
            for s in samples:
                sample2group_source.append([s, s])
        self.chart_tree("wgcna", ".sample_tree", sample_tree=newick_tree, sample2group_source = sample2group_source, json_mode="wgcna.prepare_cluster.tree.json")

    def chart_wgcna_module_tree(self, module_tree_table, group_dict=None):
        linkage_pd = pd.read_table(module_tree_table, header=None)
        linkage_list = [list(x) for x in linkage_pd.as_matrix()]
        ordered_leaf = pd.read_table(module_tree_table[:-3] + 'order.txt', header=None)[0]
        ordered_seqs = [x.strip() for x in open(module_tree_table[:-3] + 'labels.txt')]
        newick_tree = convert2(module_tree_table, ordered_leaf, ordered_seqs)
        source = list()

        module2group_source = [["name", "group"]]
        if group_dict:
            for g,ss in group_dict.items():
                for s in ss:
                    module2group_source.append([s, g])
        else:
            modules = ordered_seqs
            for s in modules:
                module2group_source.append([s, s])

        self.chart_tree("wgcna", ".module_tree", newick_tree, sample2group_source = module2group_source, json_mode="wgcna.module_cluster.tree.json")


    def chart_wgcna_module_corr(self, module_corr, module_tree):
        with open(module_tree, 'r') as tree_f:
            newick_tree = tree_f.readline().strip()[:-1]
        source = list()
        cluster_pd = pd.read_table(module_corr, index_col=0, header=0)
        ordered_seqs = list(cluster_pd.columns)
        cluster_pd["seq_id"] = cluster_pd.index
        corr_heat = [[""] + ordered_seqs]
        for line_dict in cluster_pd.iterrows():
            corr_heat.append(
                [line_dict[0]] + [line_dict[1][s] for s in ordered_seqs]
            )

        module2group_source = [["name", "group"]]
        for m in ordered_seqs:
            module2group_source.append([m, m])
        self.chart_heat_tree("wgcna", ".mofule_corr", corr_heat, newick_tree, newick_tree, module2group_source, module2group_source, "wgcna.module_relation.heatmap_tree.json")

    def denovo_chart_wgcna_module_corr(self, module_corr, module_tree):
        with open(module_tree, 'r') as tree_f:
            newick_tree = tree_f.readline().strip()[:-1]
        labels = re.findall('[(,]([^(]*?):', newick_tree)
        rows = labels
        top_group_colors = [i.lstrip('ME') for i in rows]
        heatmap_data = list()
        cluster_pd = pd.read_table(module_corr, index_col=0, header=0)
        cluster_pd_1 = cluster_pd.ix[:, labels]
        cluster_pd_2 = cluster_pd_1.reindex(index=labels)
        corr_heat = list()
        for i in labels:
            corr_heat.append(cluster_pd_2[i].tolist())

        # ordered_seqs = list(cluster_pd.columns)
        # cluster_pd["seq_id"] = cluster_pd.index
        # # corr_heat = [[""] + ordered_seqs]
        # for line_dict in cluster_pd.iterrows():
        #     corr_heat.append(
        #         # [line_dict[0]] + [line_dict[1][s] for s in ordered_seqs]
        #         [line_dict[1][s] for s in ordered_seqs]
        #     )

        self.denovo_chart_heat_tree("wgcna", ".module_corr", corr_heat, newick_tree, newick_tree, labels, top_group_colors, "wgcna_module_corr_heat.json")

    def chart_wgcna_relation_corr(self, relation_corr, relation_corr_pvalue, module_stat):
        cluster_pd = pd.read_table(relation_corr, index_col=0, header=0)
        ordered_seqs = list(cluster_pd.columns)

        visual_color = [str(colors.to_hex(c[2:]).upper()) for c in cluster_pd.index]
        cluster_pd["module"] = cluster_pd.index
        corr_heat = [[""] + ordered_seqs]
        for line_dict in cluster_pd.iterrows():
            corr_heat.append(
                [line_dict[0]] + list([float("%0.4f" %line_dict[1][x]) for x in ordered_seqs])
            )
        
        cluster_pd_c = pd.read_table(relation_corr_pvalue, index_col=0, header=0)
        ordered_seqs = list(cluster_pd_c.columns)
        cluster_pd_c["module"] = cluster_pd_c.index
        corr_p = [[""] + ordered_seqs]
        for line_dict in cluster_pd_c.iterrows():
            corr_p.append(
                [line_dict[0]] + list([float("%0.4f" %line_dict[1][x]) for x in ordered_seqs])
            )

        group_source = [["name", "group"]]
        size_stat_pd = pd.read_table(module_stat, header=0)
        for s,m in zip(size_stat_pd["size"], size_stat_pd["module"]):
            group_source.append([s, "ME" + m])

        # print corr_heat
        # print group_source
        # print corr_p
        self.chart_heat("wgcna", ".relation_heat", corr_heat, group_source, corr_p, visual_color, "wgcna.module_trait.heat.json")

    def chart_wgcna_module_column(self, module_stat):
        size_stat_pd = pd.read_table(module_stat, header=0)

        source = [
            ["item"] + list(size_stat_pd["module"]),
            ["series"] + list(size_stat_pd["size"]),
            ["category"] + list(size_stat_pd["module"])
        ]

        visual_color = [str(colors.to_hex(c).upper()) for c in size_stat_pd["module"]]
        # print source
        # print visual_color
        self.chart_column_color("wgcna", ".module_stat", source, visual_color, "wgcna.module_stat.column.json")

    def chart_wgcna_relation_ms(self, gene_trait_corr, seq_annot, module_trait_corr, module_trait_corr_pvalue):
        seq_annot_pd = pd.read_table(seq_annot,index_col=0,header=0)
        corr_pd = pd.read_table(gene_trait_corr,index_col=0,header=0)
        corr_pd.index.name = 'seq_id'
        corr_all = corr_pd.join(seq_annot_pd)
        final_corr = corr_all.reset_index()

        corr = pd.read_table(module_trait_corr, index_col=0, header=0)
        module_list = list(corr.index)
        trait_list = list(corr.columns)

        corr_p = pd.read_table(module_trait_corr_pvalue, index_col=0, header=0)
        
        for module, module_df in corr_all.groupby("module"):
            for trait in trait_list:
                source = [["name", "x", "y", "value", "category"]]
                # module_df.to_csv("test.tsv",sep = "\t", header=True)
                for seq_id, gene_dict in module_df.to_dict("index").items():
                    source.append([
                        seq_id,
                        gene_dict["kme"],
                        gene_dict[trait],
                        seq_id,
                        "all"
                    ])
                x_title = "Module Membership in {} module".format(module)
                y_title = "Gene significance"
                title = "Module membership vs. gene significance (cor={}, p<{})".format(corr.loc["ME" + module, trait], corr_p.loc["ME" + module, trait])
                self.chart_scatter("wgcna", "relation_ms." + module + "." + trait, source, x_title, y_title, title, "wgcna.module_trait_ms_gs.scatter.json")

        corr_all = corr_pd.abs().join(seq_annot_pd.loc[:, ["module"]])
        mean_gs = corr_all.groupby("module").mean()
        gs_std = corr_all.groupby("module").std()
        gs_all = mean_gs.join(gs_std, lsuffix="_gs", rsuffix="_gs_std")
        gs_all.index = ["ME"+x for x in gs_all.index]

        corr_p = pd.read_table(module_trait_corr_pvalue, index_col=0, header=0)
        corr_all2 = corr.join(corr_p, lsuffix='_corr', rsuffix="_pvalue")
        corr_all2.index.name = "module"

        module_size_pd = seq_annot_pd.groupby("module").count().iloc[:, [0]]
        module_size_pd.columns = ["size"]
        module_size_pd.index = ["ME"+x for x in module_size_pd.index]
        corr_all2 = corr_all2.join(module_size_pd)
        corr_all2 = corr_all2.join(gs_all)
        corr_all2.reset_index(inplace=True)
        corr_all2_dict_list = corr_all2.to_dict("records")
        for trait in trait_list:
            source = [
                ["item"] + list(corr_all2["module"]),
                ["series"] + list(corr_all2[trait + "_gs"]),
                ["category"] + list(corr_all2["module"])
            ]

            source2 = [
                ["item"] + list(corr_all2["module"]),
                ["series"] + [[g, s, s, ""] for g,s in zip(list(corr_all2[trait + "_gs"]), list(corr_all2[trait + "_gs_std"]))],
                ["category"] + list(corr_all2["module"])
            ]

            visual_color = [colors.to_hex(c[2:]).upper() for c in module_size_pd.index]
            # print source
            print source2
            # print visual_color

            self.chart_column_and_conf("wgcna", ".{}.relation_ms".format(trait), source, source2, visual_color, "wgcna.module_trait_ms.box.json")

    def chart_ppi_stat(self, centrality, degree):
        centrality_pd = pd.read_table(centrality,header=0)
        centrality_source = [
            list(centrality_pd["Degree_Centrality"]),
            list(centrality_pd["Closeness_Centrality"]),
            list(centrality_pd["Betweenness_Centrality"])
        ]

        categories = list(centrality_pd["Node_ID"])

        self.chart_line_highchart("ppi", ".centrality", centrality_source, categories, "ppi.center.line.json")

        degree_pd = pd.read_table(degree, header=0)
        degree_source = [list(degree_pd["Node_Num"])]
        categories = list(degree_pd["Degree"])

        self.chart_line_highchart("ppi", ".degree", degree_source, categories, "ppi.degree.line.json")

    def chart_stem_cluster(self, profile_list, cluster_table, method):
        if method == 'SCM':
            cluster = pd.read_table(cluster_table, header=0, index_col=None, sep='\t', dtype=object)

            cluster_new = cluster[['profile', 'model', 'cluster', 'pvalue']]
            cluster_list = cluster_new.to_dict('records')
            self.chart_stem("stem", ".trend", cluster_list, "stem_cluster.json")
            for i in profile_list:
                profile = pd.read_table(i, header=0, index_col=0, sep='\t')
                profile_num = profile['profile'].tolist()[0]
                profile_new = profile.drop('profile', axis=1)
                source_line = [{"data": list(profile_new.loc[j]), 'group': j} for j in profile_new.index.tolist()]
                category_line = profile_new.columns.tolist()
                source_scatter = [[h, m, profile_new[m][h], h] for h in profile_new.index.tolist() for m in profile_new.columns.tolist() ]
                source_scatter.insert(0, ['name', 'x', 'y', 'category'])
                title = "Profile #{}(0.0,-1.0,-1.0,0.0,0.0,-1.0)".format(profile_num)
                subtitle = "{} Genes Assigned;{} Genes Expected; p-value = {} (significant)".format(cluster.loc[profile_num]['assigned'], cluster.loc[profile_num]['expected'], cluster.loc[profile_num]['pvalue'])
                self.chart_stem_profile("stem", ".profile_{}".format(profile_num), source_line, category_line, source_scatter, title, subtitle,"stem_detail.json")
        elif method == 'K':
            cluster = pd.read_table(cluster_table, header=0, index_col=None, sep='\t', dtype=object)
            cluster_new = cluster[['cluster', 'model']]
            cluster_new.rename({'cluster': 'profile'}, inplace=True)
            cluster_new['cluster'] = -1
            cluster_list = cluster_new.to_dict('records')
            self.chart_stem("stem", ".trend", cluster_list, "stem_cluster.json")
            for i in profile_list:
                profile = pd.read_table(i, header=0, index_col=0, sep='\t')
                profile_num = profile['cluster'].tolist()[0]
                profile_new = profile.drop('cluster', axis=1)
                source_line = [{"data": list(profile_new.loc[j]), 'group': j} for j in profile_new.index.tolist()]
                category_line = profile_new.columns.tolist()
                source_scatter = [[h, m, profile_new[m][h], h] for h in profile_new.index.tolist() for m in profile_new.columns.tolist() ]
                source_scatter.insert(0, ['name', 'x', 'y', 'category'])
                self.chart_stem_profile("stem", ".profile_{}".format(profile_num), source_line, category_line, source_scatter, None, None, "stem_detail.json")

    def chart_tf_stat(self, column_stat, circ_gene, circ_trans):
        size_stat_pd = pd.read_table(column_stat, header=0)

        source = [
            ["item"] + list(size_stat_pd["family"]),
            ["series"] + list(size_stat_pd["transcript_num"]),
            ["category"] + list(size_stat_pd["family"])
        ]

        self.chart_column2("tf_stat", ".transcript", source, "tf.family_stat.column.json")

        source = [
            ["item"] + list(size_stat_pd["family"]),
            ["series"] + list(size_stat_pd["gene_num"]),
            ["category"] + list(size_stat_pd["family"])
        ]

        self.chart_column2("tf_stat", ".gene", source, "tf.family_stat.column.json")

        fams = list(size_stat_pd["family"])
        circ = pd.read_table(circ_trans, header=0)
        source_circ = list()
        for rec in circ.to_dict('records'):
            source_circ.append([rec["transcript_id"]] + [rec[term] for term in fams] + [rec["e_value"]])

        source_class = [[t, z] for t,z in zip(fams, list(size_stat_pd["transcript_num"]))]
        title = "Statistics of TF family"

        self.chart_circ("tf_circ", ".transcript", source_circ, source_class, title, "tf.family_stat.circ.json")


        fams = list(size_stat_pd["family"])
        circ = pd.read_table(circ_gene, header=0)
        source_circ = list()
        for rec in circ.to_dict('records'):
            source_circ.append([rec["gene_id"]] + [rec[term] for term in fams] + [rec["e_value"]])

        source_class = [[t, z] for t,z in zip(fams, list(size_stat_pd["gene_num"]))]
        title = "Statistics of TF family"

        self.chart_circ("tf_circ", ".gene", source_circ, source_class, title, "tf.family_stat.circ.json")

    def chart_json_batch(self, chart_json):
        if "filter_assemble_length" in chart_json:
            filter_assemble_length_file = chart_json['filter_assemble_length']
            self.chart_filter_assemble_length(filter_assemble_length_file)

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
    a = ChartAdvance_denovo()

    group_dict = {"SNU16": ["SNU16_4","SNU16_7","SNU16_8","SNU16_9","SNU16_1","SNU16_5","SNU16_2","SNU16_3","SNU16_6"], "H1581": ["H1581_9","H1581_7","H1581_8","H1581_3","H1581_6","H1581_1","H1581_2","H1581_4","H1581_5"]}
    subcluster_list = glob.glob("subcluster_*")

    # a.chart_wgcna_sample_tree(sys.argv[1], group_dict)
    # a.chart_wgcna_module_tree(sys.argv[1], group_dict=None)
    # a.chart_wgcna_module_corr(sys.argv[1], sys.argv[2])
    # a.chart_wgcna_module_column(sys.argv[1])
    # a.chart_wgcna_relation_corr(sys.argv[1], sys.argv[2], sys.argv[3])
    # a.chart_wgcna_relation_ms(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    #a.chart_wgcna_relation_ms

    # a.chart_ppi_stat(sys.argv[1], sys.argv[2])
    a.chart_tf_stat(sys.argv[1], sys.argv[2], sys.argv[3])
