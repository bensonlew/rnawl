# coding=utf-8
from __future__ import division
import xml.etree.ElementTree as ET
from biocluster.config import Config
import re
import collections
import json
from itertools import islice
from collections import defaultdict
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
from mbio.packages.ref_rna_v2.chart import Chart as chart_base
from mbio.packages.medical_transcriptome.chart.chart_diff_pipline import ChartDiffPipline as chart_base
reload(sys)
sys.setdefaultencoding('utf-8')
class Chart(chart_base):
    def __init__(self):
        """
        设置数据库，连接到mongod数据库，kegg_ko,kegg_gene,kegg_pathway_png三个collections2
        """

        """
        class: sample[A1, B1], group[A, B], classify: [samples, groups], level[G, T], kind[ref, new, all], category[mRNA, lncRNA]
        """
        super(Chart, self).__init__()

    def chart_genefusion_venn(self, gene_fusion_venn_file):
        source = []
        with open(gene_fusion_venn_file, "r") as r:
            r.readline()
            for line in r.readlines():
                line = line.strip().split("\t")
                source.append({"data": line[1].strip().split(","), "name": line[0]})
        self.chart_upset("genefusion_", "venn", source, "genefusion_venn.json")

    def chart_diffexp_stat(self, summary, cmp_list=None):
        with open(summary, "r") as f:
            lines = f.readlines()
            if len(lines) < 3:
                return

        summary_pd = pd.read_table(summary, header=[0, 1])
        levels = summary_pd.columns.levels
        labels = summary_pd.columns.labels
        summary_pd.columns = levels[0][labels[0]]

        categories = sorted(list(summary_pd.columns[1:-1]))
        up_list = list()
        down_list = list()
        for cmp1 in sorted(categories):
            stat = summary_pd[cmp1].value_counts()
            if "yes|up" in stat:
                up_list.append(stat["yes|up"])
            if "yes|down" in stat:
                down_list.append(stat["yes|down"])
        data1 = [up_list,down_list
        ]
        # 堆叠图
        data2 = [
            up_list,down_list

        ]
        self.chart_diffexp_stat_bar("all", ".diffexp_summary", data1, categories, "../medical_transcriptome/diffexp_stat.bar1.json")
        self.chart_diffexp_stat_bar2("all", ".diffexp_summary", data2, categories, "../medical_transcriptome/diffexp_stat.bar2.json")

    def chart_splice_diff_stat(self, event_stats, event_type, cmp_name="A_vs_B"):
        psi_data = self._parse_dic_file(file = event_type)
        diff_stats_data = self._parse_dic_file(file = event_stats)
        splice_type = ["SE", "RI", "A5SS", "A3SS", "MXE"]
        stat_type = {
            "JC": "JunctionCountOnly_event_id_set_no",
            "JCEC": "ReadsOnTargetAndJunctionCounts_event_id_set_no",
            "JC & JCEC": "JunctionCountOnly_and_ReadsOnTargetAndJunctionCounts_set_no",
            "JC | JCEC": "JunctionCountOnly_or_ReadsOnTargetAndJunctionCounts_set_no"
        }
        # print diff_stats_data
        for stat, f in stat_type.items():
            source = [["value", "name", "category"]]
            sum_a = 0
            for category in ["MXE","A5SS","RI","SE", "A3SS", "SE"]:
                sum_a += diff_stats_data[category + "_" + stat_type[stat]]
            for category in ["MXE","A5SS","RI","SE", "A3SS", "SE"]:
                num = diff_stats_data[category + "_" + stat_type[stat]]
                if sum_a == 0:
                    pct = 0
                else:
                    pct = float(num)/sum_a * 100
                source.append([
                   diff_stats_data[category + "_" + stat_type[stat]],
                    category,
                    "{}:{}%".format(category, ("%0.4f" %pct))
                ])
                title = "Summary of DES({})".format(cmp_name)

            stat_name = stat.replace(" & ", "and")
            stat_name = stat_name.replace(" | ", "or")
            self.chart_splice_diff_stat_pie(cmp_name + "_" + stat_name, ".diff_splice_stat", source, splice_type, title, "splice.in_stat.pie.json")

        stat_type = {
            "JC": "JunctionCountOnly",
            "JCEC": "ReadsOnTargetAndJunctionCounts",
        }

        for stat, f in stat_type.items():
            source = [
                {
                    "data": [psi_data[category + "_SAMPLE_1_" + f + "_exclusion"] for category in ["MXE","A5SS","RI","SE", "A3SS", "SE"]],
                    "name": "exc"
                },
                {
                    "data": [psi_data[category + "_SAMPLE_1_" + f + "_inclusion"] for category in ["MXE","A5SS","RI","SE", "A3SS", "SE"]],
                    "name": "inc"
                }
            ]
            title = "Pattern of differentially expressed AS"
            # title = "Summary of DES({})".format(cmp_name)
            self.chart_splice_diff_stat_column(cmp_name + "_" + stat, ".diff_splice_stat", source, splice_type, title, "../medical_transcriptome/splice.in_stat.column.json")

    def chart_asprofile_stat(self, asprofile_stat):
        splice_pd =  pd.read_table(asprofile_stat, header=0, sep='\t')
        categories = list(splice_pd.columns[:-1])
        source = []
        for category in categories:
            source.append({
                "category": category,
                "value": list(splice_pd[category])
            })

        samples = list(splice_pd["sample"])

        self.chart_column_stat("asprofile", ".stat", source, samples, "../medical_transcriptome/splice.all_stat.column.json")
        for rec in splice_pd.to_dict("records"):
            source = [["value", "name", "category"]]
            for category in categories:
                source.append([rec.get(category, 0),
                               category,
                               category + ":{}".format(rec.get(category))
                ])
            sample = rec["sample"]
            self.chart_pie("asprofile", "." + sample, source, "../medical_transcriptome/asprofile.stat.pie.json")

    def chart_diffexp_stat_bar(self, name, out, source, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".bar.js", 'w') as fo:
            a = json.loads(f.read())
            # a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["data"] = source
            a["categories"] = categories
            #
            # a = self.reset_margin(j=a, margin_type="bottom", word_list=categories)
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".bar.js", {"model": "highchart", "highchart_type": "showBar", "width": 600, "height": "400"}])

    def chart_diffexp_stat_bar2(self, name, out, source, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".bar2.js", 'w') as fo:
            a = json.loads(f.read())
            # a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["data"] = source
            a["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".bar2.js", {"model": "highchart", "highchart_type": "showBar", "width": 600, "height": "400"}])


    def chart_upset(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".upset.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".upset.js", {}])

    def prepare_chr_infos(self,fusion_infos):
        chr_length_file = fusion_infos["chr_length_path"]
        chr_df = pd.read_table(chr_length_file, header=None)
        chr_df.columns = ["chr", "chr_length"]
        chr_df["attributes"] = "chromosome"
        chr_up_dict = chr_df.to_dict("r")
        source_chr_length = []
        for chr in chr_up_dict:
            source_chr_length.append([chr["chr"], chr["chr_length"]])
        return source_chr_length

        # assemble_level_file = fusion_infos["assemble_level_file"]


    def chart_gene_fusion_ciros(self, sample_list, fusion_result_files,fusion_infos):
        source_chr_infos = self.prepare_chr_infos(fusion_infos)
        chr_length = fusion_infos["chr_length_path"]
        gene_pos = fusion_infos["gene_pos"]
        for sample, fusion_file in zip(sample_list, fusion_result_files):
            detail_sdf = pd.read_table(fusion_file)
            if detail_sdf.shape[0] >= 1:
                select_columns = ["#FusionName", "JunctionReadCount", "SpanningFragCount", "LeftGene", "LeftBreakpoint",
                                  "RightGene", "RightBreakpoint", "SpliceType", "FFPM"]
                select_df = detail_sdf[select_columns]
                select_df.rename({"#FusionName": "FusionName"}, axis=1, inplace=True)
                select_df.columns = map(str.lower, select_df.columns)
                select_df["fusion_unique_id"] = select_df["leftbreakpoint"] + select_df["rightbreakpoint"]
                chr_df = pd.read_table(chr_length, header=None)
                chr_df.columns = ["chr", "chr_length"]
                chr_dict = chr_df.set_index("chr").to_dict("index")
                pos_df = pd.read_table(gene_pos, header=None)
                pos_df.columns = ["gene_id", "chr", "start", "end"]
                pos_df = pos_df.set_index("gene_id")
                pos_dict = pos_df.to_dict("index")
                circos_select_columns = ["fusionname", "leftgene", "leftbreakpoint", "rightgene",
                                         "rightbreakpoint", "ffpm"]
                circos_select_df = select_df[circos_select_columns]
                circos_select_df[["left_chr", "left_local"]] = circos_select_df['leftbreakpoint'].str.split(':',expand=True).iloc[:,:2]
                circos_select_df[["right_chr", "right_local"]] = circos_select_df['rightbreakpoint'].str.split(':',expand=True).iloc[:, :2]
                circos_select_df[["left_start", "left_end"]] = circos_select_df["leftgene"].map(pos_dict).apply(pd.Series).loc[:, ["start", "end"]]
                circos_select_df[["right_start", "right_end"]] = circos_select_df["rightgene"].map(pos_dict).apply(pd.Series).loc[:, ["start", "end"]]
                circos_select_df["left_local"] = circos_select_df["left_local"].astype("int")
                circos_select_df["right_local"] = circos_select_df["right_local"].astype("int")
                if circos_select_df.shape[0] != 1:
                    if max(circos_select_df["ffpm"]) == min(circos_select_df["ffpm"]):
                        circos_select_df["width"] = 3
                    else:
                        circos_select_df["width"] = circos_select_df.apply(
                            lambda x: 4 * (x["ffpm"] - min(circos_select_df["ffpm"])) / (max(circos_select_df["ffpm"]) - min(circos_select_df["ffpm"])) + 1, axis=1)
                else:
                    circos_select_df["width"] = 3
                circos_select_line_df = circos_select_df.drop(["leftbreakpoint", "rightbreakpoint", "ffpm"], axis=1)
                source_line = []
                circos_select_line_df["LinkType"] = circos_select_line_df.apply(lambda x: "A" if x["left_chr"] == x["right_chr"] else "B", axis=1)
                circos_select_line_df["g1start"] = circos_select_line_df["left_start"]
                circos_select_line_df["g2chr"] = circos_select_line_df["right_chr"]
                circos_select_line_df["g1name"] = circos_select_line_df["leftgene"]
                circos_select_line_df["g2start"] = circos_select_line_df["right_start"]
                circos_select_line_df["LinkWidth"] = circos_select_line_df["width"]
                circos_select_line_df["fusion"] = circos_select_line_df["fusionname"]
                circos_select_line_df["g2name"] = circos_select_line_df["rightgene"]
                circos_select_line_df["g1chr"] = circos_select_line_df["left_chr"]
                circos_select_line_df["g2end"] = circos_select_line_df["right_start"]
                circos_select_line_df["g1end"] = circos_select_line_df["left_start"]
                final_select_columns = ["LinkType", "left_end", "g1start", "g2chr", "g1name", "g2start", "LinkWidth", "fusion", "g2name", "right_end", "g1chr", "g2end", "g1end"]
                f_circos_select_line_df = circos_select_line_df[final_select_columns]
                circos_select_line_detail = f_circos_select_line_df.to_dict("r")
                source_line = circos_select_line_detail
                source_tree = []
                all_genes = circos_select_df["leftbreakpoint"].tolist() + circos_select_df["rightbreakpoint"].tolist()
                tree_dict = defaultdict(list)
                for index, item in circos_select_df.iterrows():
                    left_gene_name = item["fusionname"].split("-")[0]
                    left_gene_chr = item["left_chr"]
                    left_gene_loc = item["left_local"]
                    left_area_loc = int(math.ceil(left_gene_loc/chr_dict[left_gene_chr]["chr_length"] * 100))
                    left_detail_info = [(left_area_loc - 1) * chr_dict[left_gene_chr]["chr_length"] / 100, left_gene_loc,
                                        left_gene_name]
                    tree_dict[left_gene_chr + "_" + str(left_area_loc)].append(left_detail_info)
                    right_gene_name = item["fusionname"].split("-")[-1]
                    right_gene_chr = item["right_chr"]
                    right_gene_loc = item["right_local"]
                    right_area_loc = int(math.ceil(right_gene_loc / chr_dict[right_gene_chr]["chr_length"] * 100))
                    right_detail_info = [(right_area_loc - 1) * chr_dict[right_gene_chr]["chr_length"] / 100,
                                         right_gene_loc, right_gene_name]
                    tree_dict[right_gene_chr + "_" + str(right_gene_loc)].append(right_detail_info)
                for cluster in tree_dict:
                    # sample_name = sample_name
                    chr = cluster.split("_")[0]
                    root = tree_dict[cluster][0][0]
                    branch = [i[1] for i in tree_dict[cluster]]
                    gene_name = [i[2] for i in tree_dict[cluster]]
                    data = {"chr":chr,
                            "root":root,
                            "genename":gene_name,
                            "branch":branch}
                    source_tree.append(data)
                categories = "categories"
                self.chart_gene_fusion(sample, "gene_fusion", source_chr_infos, source_line,source_tree,categories,
                                       "genefusion.json")

    def chart_gene_fusion(self,name, out, source_chr_infos,source_line, source_tree,categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source_chr_infos
            a["dataset"][1]["source"] = source_line
            a["dataset"][2]["source"] = source_tree
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".js", {}])

    def chart_snp_dis(self, snp_position_distribution, samples=None):
        snp_pos_pd =  pd.read_table(snp_position_distribution, header=0, sep='\t')
        snp_pos_columns = list(snp_pos_pd.columns)
        snp_pos_columns.remove("range_key")
        snp_pos_columns.remove("stat_type")
        samples = snp_pos_columns
        # samples = snp_pos_pd.columns[: -2]

        for sample in samples:
            source =  [["name","value","category"]]
            total = sum(snp_pos_pd[sample])
            for n,c in zip(snp_pos_pd[sample], snp_pos_pd["range_key"]):
                source.append([
                    n,
                    c,
                    "{}:{}%".format(c, "%0.2f" %(float(n)/total * 100))
                ])

            title = "SNP distribution in genome regions ({sample_name})".format(sample_name=sample)
            self.chart_pie(sample, ".snp.pos_stat", source, "snp.pos_stat.pie.json", title)

    def chart_somatic_dis(self, snp_position_distribution, samples=None):
        snp_pos_pd =  pd.read_table(snp_position_distribution, header=0, sep='\t')
        snp_pos_columns = list(snp_pos_pd.columns)
        snp_pos_columns.remove("range_key")
        snp_pos_columns.remove("stat_type")
        samples = snp_pos_columns
        # samples = snp_pos_pd.columns[: -2]

        for sample in samples:
            source =  [["name","value","category"]]
            total = sum(snp_pos_pd[sample])
            for n,c in zip(snp_pos_pd[sample], snp_pos_pd["range_key"]):
                source.append([
                    n,
                    c,
                    "{}:{}%".format(c, "%0.2f" %(float(n)/total * 100))
                ])

            title = "Somatic distribution in genome regions ({sample_name})".format(sample_name=sample)
            self.chart_pie(sample, ".somatic.pos_stat", source, "snp.pos_stat.pie.json", title)

    def chart_snp_stat(self, snp_stat, samples=None):
        snp_stat_pd =  pd.read_table(snp_stat, header=0, sep='\t')
        stat_pd_columns = list(snp_stat_pd.columns)
        stat_pd_columns.remove("range_key")
        stat_pd_columns.remove("stat_type")
        samples = stat_pd_columns
        for sample in samples:
            source =  [["name","value","category"]]
            total = sum(snp_stat_pd[sample])
            for n,c in zip(snp_stat_pd[sample], snp_stat_pd["range_key"]):
                source.append([
                    n,
                    c,
                    "{}:{}%".format(c, "%0.2f" %(float(n)/total * 100))
                ])

            title = "Summary of SNP ({sample_name})".format(sample_name=sample)
            self.chart_pie(sample, ".snp.type_stat", source, "snp.type_stat.pie.json", title)

            source = [{
                "data": list(snp_stat_pd[sample]),
                "name": sample
            }]

            categories = list(snp_stat_pd["range_key"])
            title = "Summary of SNP"
            self.chart_column(sample, ".snp.type_stat", source, categories, title, "snp.type_stat.column.json")

    def chart_snp_depth_stat(self, snp_depth_stat, samples=None):
        snp_stat_pd =  pd.read_table(snp_depth_stat, header=0, sep='\t')
        stat_pd_columns = list(snp_stat_pd.columns)
        stat_pd_columns.remove("range_key")
        stat_pd_columns.remove("stat_type")
        samples = stat_pd_columns
        for sample in samples:
            source =  [["name","value","category"]]
            total = sum(snp_stat_pd[sample])
            for n,c in zip(snp_stat_pd[sample], snp_stat_pd["range_key"]):
                source.append([
                    n,
                    c,
                    "{}:{}%".format(c, "%0.2f" %(float(n)/total * 100))
                ])

            title = "Summary of SNP ({sample_name})".format(sample_name=sample)
            self.chart_pie(sample, ".snp.depth_stat", source, "snp.depth_stat.pie.json", title)

            source = [{
                "data": list(snp_stat_pd[sample]),
                "name": sample
            }]

            categories = list(snp_stat_pd["range_key"])
            title = "Summary of SNP"
            self.chart_column(sample, ".snp.depth_stat", source, categories, title, "snp.depth_stat.column.json")

    def chart_somatic_stat(self, snp_stat, samples=None):
        snp_stat_pd =  pd.read_table(snp_stat, header=0, sep='\t')
        stat_pd_columns = list(snp_stat_pd.columns)
        stat_pd_columns.remove("range_key")
        stat_pd_columns.remove("stat_type")
        samples = stat_pd_columns
        for sample in samples:
            source =  [["name","value","category"]]
            total = sum(snp_stat_pd[sample])
            for n,c in zip(snp_stat_pd[sample], snp_stat_pd["range_key"]):
                source.append([
                    n,
                    c,
                    "{}:{}%".format(c, "%0.2f" %(float(n)/total * 100))
                ])

            title = "Summary of Somatic ({sample_name})".format(sample_name=sample)
            self.chart_pie(sample, ".somatic.type_stat", source, "snp.type_stat.pie.json", title)

            source = [{
                "data": list(snp_stat_pd[sample]),
                "name": sample
            }]

            categories = list(snp_stat_pd["range_key"])
            title = "Summary of Somatic"
            self.chart_column(sample, ".somatic.type_stat", source, categories, title, "snp.type_stat.column.json")

    def chart_somatic_depth_stat(self, snp_depth_stat, samples=None):
        snp_stat_pd =  pd.read_table(snp_depth_stat, header=0, sep='\t')
        # samples = snp_stat_pd.columns[: -2]
        stat_pd_columns = list(snp_stat_pd.columns)
        stat_pd_columns.remove("range_key")
        stat_pd_columns.remove("stat_type")
        samples = stat_pd_columns
        for sample in samples:
            source =  [["name","value","category"]]
            total = sum(snp_stat_pd[sample])
            for n,c in zip(snp_stat_pd[sample], snp_stat_pd["range_key"]):
                source.append([
                    n,
                    c,
                    "{}:{}%".format(c, "%0.2f" %(float(n)/total * 100))
                ])

            title = "Summary of Somatic ({sample_name})".format(sample_name=sample)
            self.chart_pie(sample, ".somatic.depth_stat", source, "snp.depth_stat.pie.json", title)

            source = [{
                "data": list(snp_stat_pd[sample]),
                "name": sample
            }]

            categories = list(snp_stat_pd["range_key"])
            title = "Summary of Somatic"
            self.chart_column(sample, ".somatic.depth_stat", source, categories, title, "snp.depth_stat.column.json")


    def chart_diffexp_scatter(self, diff_exp, cmps, soft=None,pvalue_pajust = "pajust"):
        diff_pd = pd.read_table(diff_exp, header=0, sep='\t')
        if soft.lower() == "noiseq":
            need_cols = ['seq_id', 'fc', 'log2fc', 'D', 'prob', 'significant', 'regulate']
            volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', 'D', 'significant', 'regulate']]
            volcano_pd.columns = ['seq_id', 'log2fc', 'D', 'significant', 'regulate']
            volcano_pd_nosig = volcano_pd[volcano_pd['significant'] == 'no']
            volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']

            #add by fwy 20210304
            fname = os.path.basename(diff_exp)
            ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(soft.lower()), fname).groups()
            diff_pd = pd.read_table(diff_exp, header=0, sep='\t')
            columns = diff_pd.columns
            fc_ind = list(columns).index('fc')
            need_cols = ['seq_id', '{}_mean'.format(ctrl), '{}_mean'.format(test), 'fc', 'log2fc', 'D', 'prob',
                         'significant', 'regulate']
            need_cols += [columns[fc_ind - 4], columns[fc_ind - 3]]
            print(need_cols)
            samples = list()
            for x in columns:
                _m = re.match(r'(.*)_count$', x)
                if _m:
                    samples.append(_m.groups()[0])

            cmp_combine = ctrl + '|' + test
            cmp_pd = pd.DataFrame([cmp_combine] * diff_pd.shape[0], columns=['compare'])
            tmp_pd = pd.concat([diff_pd.loc[:, need_cols], cmp_pd], axis=1)
            tmp_pd.columns = list(tmp_pd.columns[:-3]) + ['group1', 'group2', 'compare']


            if volcano_pd_nosig.shape[0] > 8000:
                volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
                volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.5)
                if volcano_pd_nosig.shape[0] > 12000:
                    volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
                if volcano_pd_nosig.shape[0] > 12000:
                    volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
                    volcano_pd = pd.concat([volcano_pd_sig, volcano_pd_nosig], axis=0)
            volcano_pd_nosig.rename(columns={"D": "log10pvalue"}, inplace=True)
            volcano_pd_sig.rename(columns={"D": "log10pvalue"}, inplace=True)
        else:
            columns = diff_pd.columns
            fc_ind = list(columns).index('fc')
            need_cols = ['seq_id', 'fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
            need_cols += [columns[fc_ind - 2], columns[fc_ind - 1]]
            print(need_cols)
            samples = list()
            for x in columns:
                _m = re.match(r'(.*)_count$', x)
                if _m:
                    samples.append(_m.groups()[0])
            fname = os.path.basename(diff_exp)
            ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(soft.lower()), fname).groups()
            cmp_combine = ctrl + '|' + test
            cmp_pd = pd.DataFrame([cmp_combine] * diff_pd.shape[0], columns=['compare'])
            tmp_pd = pd.concat([diff_pd.loc[:, need_cols], cmp_pd], axis=1)
            tmp_pd.columns = list(tmp_pd.columns[:-3]) + ['group1', 'group2', 'compare']

            pvalue_padjust = pvalue_pajust
            status_list, stat_cutoff = self._get_volcano_status_cutoff(diff_pd, pvalue_padjust)
            #print "stat_cutoff is {}".format(stat_cutoff)
            volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', pvalue_padjust, 'significant', 'regulate']]
            bool_ind = volcano_pd[pvalue_padjust] <= 0
            min_pvalue = min([x if x > 0 else '' for x in volcano_pd[pvalue_padjust].tolist()])
            volcano_pd.loc[bool_ind, pvalue_padjust] = min_pvalue
            volcano_pd[pvalue_padjust] = -np.log10(volcano_pd[pvalue_padjust])
            volcano_pd.dropna(inplace=True)
            volcano_pd.columns = ['seq_id', 'log2fc', 'log10pvalue', 'significant', 'regulate']
            bool_ind = volcano_pd['log10pvalue'] > stat_cutoff
            # volcano_pd.loc[bool_ind, 'log10pvalue'] = stat_cutoff
            volcano_pd_nosig = volcano_pd[volcano_pd['significant'] == 'no']
            volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
            if volcano_pd_nosig.shape[0] > 8000:
                volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
                volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.5)
                if volcano_pd_nosig.shape[0] > 12000:
                    volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
                if volcano_pd_nosig.shape[0] > 12000:
                    volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
                volcano_pd = pd.concat([volcano_pd_sig, volcano_pd_nosig], axis=0)

        down_list = list()
        up_list = list()
        nosig_list = list()

        down_list_s = list()
        up_list_s = list()
        nosig_list_s = list()
        for vol in volcano_pd_sig.to_dict('records'):
            if vol["regulate"] == "up":
                up_list.append({
                    "y": vol["log10pvalue"],
                    "x": float(vol["log2fc"]),
                    "selected": False,
                    "name": vol["seq_id"],
                    "color": "#FF2020"
                })
            else:
                down_list.append({
                    "y": vol["log10pvalue"],
                    "x": float(vol["log2fc"]),
                    "selected": False,
                    "name": vol["seq_id"],
                    "color": "#388E3C"
                })
        for vol in volcano_pd_nosig.to_dict('records'):
            nosig_list.append({
                "y": vol["log10pvalue"],
                "x": float(vol["log2fc"]),
                "selected": False,
                "name": vol["seq_id"],
                "color": "#808080"
            })

        data = {
            "down": down_list,
            "up": up_list,
            "nosig": nosig_list
        }

        title = "{}.volcano".format(cmps)
        y_title = "-Log10(Padjust)"
        if pvalue_pajust.lower() == "pvalue" :
            y_title = "-Log10(Pvalue)"
        if soft.lower() == "noiseq":
            self.chart_diffexp_volcano(cmps, ".diffexp", data, title, "diffexp_volcano.scatter.json", x_title="Log2FC", y_title="D")
        else:
            self.chart_diffexp_volcano(cmps, ".diffexp", data, title, "diffexp_volcano.scatter.json",y_title = y_title )

        # add ma data 20200824
        ma_pd = tmp_pd.loc[:, ['seq_id', 'group1', 'group2', 'compare', 'significant', 'regulate', 'log2fc']]
        ma_pd.set_index('seq_id', inplace=True)
        ma_pd = ma_pd[ma_pd["regulate"] != "no test"]
        ma_pd = ma_pd.loc[volcano_pd['seq_id'], :].reset_index()
        ma_pd['mean_tpm'] = (ma_pd['group1'] + ma_pd['group2'] + 1).apply(np.log10)
        ma_pd = ma_pd.drop(columns=['group1', 'group2'])
        # ma_pd['group1'] = (ma_pd['group1'] + 1).apply(np.log10)
        # ma_pd['group2'] = (ma_pd['group2'] + 1).apply(np.log10)
        # ma_dict_list += ma_pd.to_dict('records')


        for vol in ma_pd.to_dict('records'):
            if vol["regulate"] == "up" and vol["significant"] == "yes":
                up_list_s.append({
                    "y": vol["log2fc"],
                    "x": vol["mean_tpm"],
                    "selected": False,
                    "name": vol["seq_id"],
                    "color": "#FF2020"
                })
            elif vol["regulate"] == "down" and vol["significant"] == "yes":
                down_list_s.append({
                    "y": vol["log2fc"],
                    "x": vol["mean_tpm"],
                    "selected": False,
                    "name": vol["seq_id"],
                    "color": "#388E3C"
                })
            else:
                nosig_list_s.append({
                    "y": vol["log2fc"],
                    "x": vol["mean_tpm"],
                    "selected": False,
                    "name": vol["seq_id"],
                    "color": "#808080"
                })

        data = {
            "down": down_list_s,
            "up": up_list_s,
            "nosig": nosig_list_s
        }
        def get_exp_type(df):
            for col in df.columns:
                if col.lower().endswith("tpm"):
                    return "tpm"
            else:
                return "fpkm"
        exp_tpye = get_exp_type(diff_pd)
        title = "{}.MA plot".format(cmps)
        x_label = "log10({})".format(exp_tpye)
        y_label = "log2(FC)"
        self.chart_diffexp_scatter3(cmps, ".diffexp_ma", data, title, x_label, y_label, "../medical_transcriptome/diffexp_ma.scatter.json")

    def chart_diffexp_scatter3(self, name, out, data, title, x_label, y_label, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".scatter.js", 'w') as fo:
            a = json.loads(f.read())
            a["params"]["x_label"] = x_label
            a["params"]["y_label"] = y_label
            a["params"]["title"] = title
            a["data"] = data
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".scatter.js", {"model": "highchart","use_medical":"yes"}])

    def chart_diff_genesets_venn(self, diff_genesets_outdir,genesets):
        # 通过获取差异基因集结果文件的原始ids文件来获取信息
        source_venn = list()
        for geneset in genesets:
            geneset_dir = os.path.join(diff_genesets_outdir,geneset)
            geneset_info = json.load(open(os.path.join(geneset_dir, "analysis_json")))
            #modify by fwy 20210715 0个基因的基因集不在绘制venn图
            if int(geneset_info["gene_num"]) >0:
                genset_list_file = geneset_info["geneset_list"]
                with open(genset_list_file,"r") as r:
                    genes = r.read().strip().split("\n")
                source_venn.append({
                    "data": genes,
                    "name": geneset
                })
        self.chart_venn("diff_genesets", ".analysis", source_venn, "geneset.venn.venn.json")
        if len(genesets) <= 6:
            self.chart_upset("diff_genesets", ".analysis", source_venn, "geneset.venn.upset.json")
        # self.chart_exp_venn_venn("diff_genesets", ".analysis", source_venn, "exp.relation.venn.json")
        # self.chart_exp_upset("diff_genesets", ".analysis", source_venn, "annotation.stat.upset.json")

    def get_seqid2name(self,seqid2namefile):
        annot_info = pd.read_table(seqid2namefile)
        gene2name_df = annot_info[["transcript_id", "gene_id", "gene_name"]]
        gene2name_df["gene_name"].fillna(gene2name_df["gene_id"], inplace=True)
        id2namedict = {}
        id2namedict = dict(zip(gene2name_df['gene_id'], [x if x else '-' for x in gene2name_df['gene_name']]))
        return id2namedict


    def chart_geneset_cluster(self, cluster_exp, cluster_tree, sample_tree, subcluster_list, group_dict=None,
                              samples_order=None,seqid2namefile = None):
        seq_id2name = self.get_seqid2name(seqid2namefile)
        with open(sample_tree, 'r') as f:
            sample_tree = f.readline().strip()
            samples = f.readline().strip().split(";")

        with open(cluster_tree, 'r') as f:
            gene_tree = f.readline().strip()
            genes = f.readline().strip().split(";")

        cluster_pd = pd.read_table(cluster_exp, index_col=0, header=0)
        gene_num = cluster_pd.shape[0]
        chart_height = 600 if gene_num <= 600 else 1500

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
            for g, ss in group_dict.items():
                for s in ss:
                    # sample2group_dict
                    sample2group_source.append([s, g])
        else:
            for s in samples:
                sample2group_source.append([s, s])

        gene2group_source = [["name", "group"]]
        gene2group_dict = dict()
        for c in subcluster_list:
            print
            c
            c_num = c.split("subcluster_")[1].split(".")[0].split("_")[0]
            sub_genes = list()
            with open(c, 'r') as f:
                f.readline()
                for line in f:
                    sub_genes.append(line.strip().split()[0])
                    gene2group_dict[line.strip().split()[0]] = int(c_num)
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
            self.chart_line_point("cluster", "." + c_num, source_line, source_point, categories, title,
                                  "geneset.cluster.line.json")

        for gene in genes:
            gene2group_source.append([gene, gene2group_dict[gene]])

        self.chart_heat_tree("geneset", ".cluster", corr_heat, sample_tree, gene_tree, sample2group_source,
                             gene2group_source, "geneset.cluster.heatmap.json",height=chart_height)

    def chart_exp_venn(self, exp_venn_file):
        # 使用参考转录本做
        with open(exp_venn_file, "r") as fr:
            source_venn = list()
            fr.readline()
            for line in fr:
                cols = line.strip().split("\t")
                # print cols
                source_venn.append({
                    "data": cols[1].split(","),
                    "name": cols[0]
                })

        self.chart_exp_venn_venn("all", ".exp", source_venn, "exp.relation.venn.json")
        self.chart_exp_upset("all", ".exp", source_venn, "annotation.stat.upset.json")


    def chart_json_batch(self, chart_json):
        if "qc_file_raw" in chart_json:
            qc_files = [chart_json["qc_file_raw"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.chart_raw_qc(chart_json["samples"], qc_files)

        if "qc_file_use" in chart_json:
            qc_files = [chart_json["qc_file_use"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.chart_raw_qc(chart_json["samples"], qc_files, qc_type="clean")

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

        if "assemble_stat" in chart_json:
            ref_gene2tran, ref_tran2exon, new_gene2tran, new_tran2exon = None, None, None, None

            if os.path.exists(chart_json["assemble_stat"] + "/old_genes.gtf.trans_1.txt"):
                ref_gene2tran = chart_json["assemble_stat"] + "/old_genes.gtf.trans_1.txt"

            if os.path.exists(chart_json["assemble_stat"] + "/old_trans.gtf.exon_1.txt"):
                ref_tran2exon = chart_json["assemble_stat"] + "/old_trans.gtf.exon_1.txt"

            if os.path.exists(chart_json["assemble_stat"] + "/new_genes.gtf.trans_1.txt"):
                new_gene2tran = chart_json["assemble_stat"] + "/new_genes.gtf.trans_1.txt"

            if os.path.exists(chart_json["assemble_stat"] + "/new_trans.gtf.exon_1.txt"):
                new_tran2exon = chart_json["assemble_stat"] + "/new_trans.gtf.exon_1.txt"

            self.chart_assemble_relation(ref_gene2tran, ref_tran2exon, new_gene2tran, new_tran2exon)

        if "annot_stat" in chart_json:
            annot_stat = chart_json["annot_stat"]
            gene_exp = chart_json["gene_exp_all"]
            trans_exp = chart_json["tran_exp_all"]
            venn_dir = os.path.dirname(chart_json["annot_stat"])
            self.chart_annotation_stat("", gene_exp, trans_exp, annot_stat, venn_dir)

        if "gene_exp_ref" in chart_json:
            gene_exp = chart_json["gene_exp_ref"]
            tran_exp = chart_json["tran_exp_ref"]
            group_dict = chart_json["group_dict"]
            self.chart_exp_dis("", gene_exp, tran_exp, group_dict)

        if "exp_venn" in chart_json and chart_json["exp_venn"]:
            venn = chart_json["exp_venn"]
            self.chart_exp_venn(venn)

        if "exp_pca_file" in chart_json:
            group_dict = chart_json["group_dict"]
            exp_pca_file = chart_json["exp_pca_file"]
            exp_pca_var_file = chart_json["exp_pca_var_file"]
            if "exp_pca_ellipse" in chart_json:
                exp_pca_ellipse = chart_json["exp_pca_ellipse"]
            else:
                exp_pca_ellipse = None
            self.chart_exp_pca(exp_pca_file, exp_pca_var_file, group_dict=group_dict, exp_pca_ellipse=exp_pca_ellipse, pcs=["PC1", "PC2"])

        if "exp_corr_file" in chart_json:
            group_dict = chart_json["group_dict"]
            exp_corr_file = chart_json["exp_corr_file"]
            exp_corr_tree_file = chart_json["exp_corr_tree_file"]
            # exp_pca_ellipse = sys.argv[3]
            self.chart_exp_corr(exp_corr_tree_file, exp_corr_file, group_dict)

        #同有参,医学差异模块图片在tool中生成
        '''
        if "diff_exp_summary" in chart_json:
            diff_exp_summary = chart_json["diff_exp_summary"]
            self.chart_diffexp_stat(diff_exp_summary)

        if "diff_exp" in chart_json:
            for cmps in chart_json["cmp_list"]:
                diff_exp = chart_json["diff_exp"].format(control=cmps[0], test=cmps[1])
            self.chart_diffexp_scatter(diff_exp, "{}_vs_{}".format(cmps[0], cmps[1]))
        '''

        if "splice_stat" in chart_json:
            for splicestat in ["JC", "JCEC"]:
                splice_stat = chart_json["splice_stat"].format(splicestat=splicestat)
                self.chart_splice_all_stat(splice_stat,  splicestat)

        if "splice_diff" in chart_json:
            for cmps in chart_json["cmp_list"]:
                splice_diff =  chart_json["splice_diff"].format(control=cmps[0], test=cmps[1])
                splice_psi = chart_json["splice_psi"].format(control=cmps[0], test=cmps[1])
                self.chart_splice_diff_stat(splice_diff, splice_psi, cmp_name="_vs_".join(cmps))

        if "as_profile_stat" in chart_json:
            self.chart_asprofile_stat(chart_json["as_profile_stat"])

        if "snp_distribution" in chart_json:
            self.chart_snp_dis(chart_json["snp_distribution"])
        if "snp_stat" in chart_json:
            self.chart_snp_stat(chart_json["snp_stat"])
        if "snp_depth" in chart_json:
            self.chart_snp_depth_stat(chart_json["snp_depth"])

        #医学专属somatic和snp区别仅在于文件名
        if "somatic_distribution" in chart_json:
            self.chart_somatic_dis(chart_json["somatic_distribution"])
        if "somatic_stat" in chart_json:
            self.chart_somatic_stat(chart_json["somatic_stat"])
        if "somatic_depth" in chart_json:
            self.chart_somatic_depth_stat(chart_json["somatic_depth"])

        if "gene_fusion_venn" in chart_json:
            self.chart_genefusion_venn(chart_json["gene_fusion_venn"])

        if "genefusion" in chart_json:
            fusion_files = [chart_json["genefusion"]["result_dir"].format(sample_name = sample) for sample in chart_json["samples"]]
            self.chart_gene_fusion_ciros(chart_json["samples"],fusion_files,chart_json["genefusion"])

        #以下为差异一键化工作流
        if  "diff_geneset_venn" in chart_json:
            try:
                diff_genesets_out =  chart_json["diff_geneset_venn"]
                self.chart_diff_genesets_venn(diff_genesets_out,chart_json["genesets"])
            except:
                pass

        if "cluster_geneset_name" in chart_json:
            self.chart_geneset_cluster(chart_json["cluster_exp"], chart_json["cluster_tree"], chart_json["sample_tree"],
                                       chart_json["subcluster_list"], group_dict=chart_json["group_dict"],seqid2namefile = chart_json["gene_annot_file"] )


        if "go_class" in chart_json:
            go_class_files = [chart_json["go_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_go_class(chart_json["genesets"], go_class_files)

        if "go_enrich" in chart_json:
            go_enrich_files = [chart_json["go_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_go_enrich(chart_json["genesets"], go_enrich_files)

        if "kegg_class" in chart_json:
            kegg_class_files = [chart_json["kegg_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_kegg_class(chart_json["genesets"], kegg_class_files,chart_json["kegg_level"])

        if "kegg_enrich" in chart_json:
            kegg_enrich_files = [chart_json["kegg_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_kegg_enrich(chart_json["genesets"], kegg_enrich_files)

        if "reactome_class" in chart_json:
            reactome_class_files = [chart_json["reactome_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_reactome_class(chart_json["genesets"], reactome_class_files)

        if "reactome_enrich" in chart_json:
            reactome_enrich_files = [chart_json["reactome_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_reactome_enrich(chart_json["genesets"], reactome_enrich_files)

        if "do_class" in chart_json:
            do_class_files = [chart_json["do_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_do_class(chart_json["genesets"], do_class_files)

        if "do_enrich" in chart_json:
            do_enrich_files = [chart_json["do_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_do_enrich(chart_json["genesets"], do_enrich_files)

        self.generate_html_sh()

if __name__ == '__main__':
    a = Chart()


    chart_json = {
        "samples": ["shNS_1", "shNS_2", "shNS_3", "shF2_1", "shF2_2", "shF2_3", "shA1_1", "shA1_2", "shA1_3", "DMSO_1", "DMSO_2", "DMSO_3","YC49_1","YC49_2","YC49_3","ZC01_1","ZC01_2","ZC01_3"],
        "group_dict": {"shNS_87": ["shNS_87", "shNS_87", "shNS_87"], "shF2_87": ["shF2_1", "shF2_2", "shF2_3"], "shA1_87": ["shA1_1", "shA1_2", "shA1_3"],
                       "DMSO_87": ["DMSO_1", "DMSO_2", "DMSO_3"],"YC49_87":["YC49_1","YC49_2","YC49_3"],"ZC01_87":["ZC01_1","ZC01_2","ZC01_3"]},
        "genefusion":{"result_dir":"/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test_20210202/genefusion/output/star_fusion/{sample_name}/star-fusion.fusion_predictions.abridged.tsv",
                        "chr_length_path":"/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test_20210202/genefusion/chr_length",
                        "gene_pos":"/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test_20210202/genefusion/gene_pos"}}

    a.chart_json_batch(chart_json)

    chart_json = {
        "samples": ["A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3"],
        "group_dict": {"A":["A1","A2","A3"], "B": ["B1","B2","B3"], "C": ["C1","C2","C3"], "D": ["D1","D2","D3"]},
        "qc_file": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/HiseqReadsStat/output/qualityStat/{sample_name}.l.qual_stat,/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/HiseqReadsStat/output/qualityStat/{sample_name}.r.qual_stat",
        "align_satu_r": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/saturation/satur_{sample_name}.eRPKM.xls.saturation.R",
        "align_satu_p": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/saturation/satur_{sample_name}.eRPKM.xls.cluster_percent.xls",
        "align_coverage": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/coverage/{sample_name}.geneBodyCoverage.txt",
        "align_pos": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/distribution/{sample_name}.reads_distribution.txt",
        "align_chr": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/chr_stat/{sample_name}.bam_chr_stat.xls",
        "assemble_step": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/RefrnaAssemble/RefassembleStat/output/trans_count_stat_200.txt",
        "assemble_new": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/RefrnaAssemble/RefassembleStat/output/code_num.txt",
        "gene_exp_ref": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/ref.gene.tpm.matrix",
        "tran_exp_ref": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/ref.transcript.fpkm.matrix",
        "gene_exp_all": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/gene.tpm.matrix",
        "tran_exp_all": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/transcript.fpkm.matrix",
        "annot_stat": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/AnnotMerge__1/output/allannot_class/all_stat.xls"}



    # chart_json = {
    #     "samples": ["A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3"],
    #     "group_dict": {"A":["A1","A2","A3"], "B": ["B1","B2","B3"], "C": ["C1","C2","C3"], "D": ["D1","D2","D3"]},
    #     "qc_file": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/HiseqReadsStat/output/qualityStat/{sample_name}.l.qual_stat,/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/HiseqReadsStat/output/qualityStat/{sample_name}.r.qual_stat",
    #     "align_satu_r": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/saturation/satur_{sample_name}.eRPKM.xls.saturation.R",
    #     "align_satu_p": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/saturation/satur_{sample_name}.eRPKM.xls.cluster_percent.xls",
    #     "align_coverage": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/coverage/{sample_name}.geneBodyCoverage.txt",
    #     "align_pos": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/distribution/{sample_name}.reads_distribution.txt",
    #     "align_chr": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/chr_stat/{sample_name}.bam_chr_stat.xls",
    #     "assemble_step": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/RefrnaAssemble/RefassembleStat/output/trans_count_stat_200.txt",
    #     "assemble_new": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/RefrnaAssemble/RefassembleStat/output/code_num.txt",
    #     "gene_exp_ref": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/ref.gene.tpm.matrix",
    #     "tran_exp_ref": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/ref.transcript.fpkm.matrix",
    #     "gene_exp_all": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/gene.tpm.matrix",
    #     "tran_exp_all": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/transcript.fpkm.matrix",
    #     "annot_stat": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/AnnotMerge__1/output/allannot_class/all_stat.xls"}
    #

    # a.chart_json_batch(chart_json)
    # a.chart_satu(samples, satu_r_file, satu_p_file)

    # cov_files = sys.argv[2].split(";")
    # a.chart_coverage(samples, cov_files)

    # pos_files = sys.argv[2].split(";")
    # a.chart_readpos(samples, pos_files)

    # assemble_new_files = sys.argv[2]
    # a.chart_assemble_new(assemble_new_files)

    # assemble_step_files = sys.argv[2]
    # a.chart_assemble(assemble_step_files)

    # assemble_step_files = sys.argv[2]
    # a. (assemble_step_files)

    # names = sys.argv[1].split(";")
    # annot_stat = sys.argv[2]
    # gene_exp = sys.argv[3]
    # trans_exp = sys.argv[4]
    # venn_dir = sys.argv[5]
    # a.chart_annotation_stat(names, gene_exp, trans_exp, annot_stat, venn_dir)

    # gene_exp = sys.argv[2]
    # tran_exp = sys.argv[3]
    # group_dict = {"A":["A1","A2","A3"], "B": ["B1","B2","B3"], "C": ["C1","C2","C3"], "D": ["D1","D2","D3"]}
    # a.chart_exp_dis(samples, gene_exp, tran_exp, group_dict)


    # venn = sys.argv[1]
    # a.chart_exp_venn(venn)

    # group_dict = {"A":["A1","A2","A3"], "B": ["B1","B2","B3"], "C": ["C1","C2","C3"], "D": ["D1","D2","D3"]}
    # exp_pca_file = sys.argv[1]
    # exp_pca_var_file = sys.argv[2]
    # exp_pca_ellipse = sys.argv[3]
    # a.chart_exp_pca(exp_pca_file, exp_pca_var_file, group_dict=group_dict, exp_pca_ellipse=exp_pca_ellipse, pcs=["PC1", "PC2"])

    # group_dict = {"A":["A1","A2","A3"], "B": ["B1","B2","B3"], "C": ["C1","C2","C3"], "D": ["D1","D2","D3"]}
    # exp_corr_file = sys.argv[2]
    # exp_corr_tree_file = sys.argv[1]
    # # exp_pca_ellipse = sys.argv[3]
    # a.chart_exp_corr(exp_corr_tree_file, exp_corr_file, group_dict)

    # diff_exp_summary = sys.argv[1]
    # a.chart_diffexp_stat(diff_exp_summary)

    # diff_exp = sys.argv[1]
    # a.chart_diffexp_scatter(diff_exp, "A_vs_B")

    # splice_stat = sys.argv[1]
    # a.chart_splice_all_stat(splice_stat)

    # splice_diff =  sys.argv[1]
    # splice_psi = sys.argv[2]
    # a.chart_splice_diff_stat(splice_diff, splice_psi)

    # snp_dis= sys.argv[1]
    # snp_stat = sys.argv[2]
    # a.chart_snp_dis(snp_dis)
    # a.chart_snp_stat(snp_stat)

    # a.chart_assemble_relation(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    # a.chart_asprofile_stat(sys.argv[1])
