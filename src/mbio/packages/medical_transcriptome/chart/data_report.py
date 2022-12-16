# coding=utf-8
from __future__ import division
from curses import raw
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
class DataReport(chart_base):
    def __init__(self):
        """
        设置数据库，连接到mongod数据库，kegg_ko,kegg_gene,kegg_pathway_png三个collections2
        """

        """
        class: sample[A1, B1], group[A, B], classify: [samples, groups], level[G, T], kind[ref, new, all], category[mRNA, lncRNA]
        """
        self.png_json = {}
        super(DataReport, self).__init__()



    def data_frame_to_json(self, df, rename=dict()):
        columns = df.columns
        table_header = [{"label":x, "field":x} for x in columns]
        table = []
        for row in df.iterrows():
            rowdict = dict(row[1])
            table.append(rowdict)
        return table_header, table

    def run(self, data_json, png_json=None):
        with open(data_json, 'r') as f:
            a = json.load(f)
        report_json =self.chart_json_batch(a)

        if png_json:
            with open(png_json, 'r') as f:
                b = json.load(f)
        report_png = self.png_dict(b)

        report_json["png"] = report_png

        with open("report.js", "w") as f:
            json_str = json.dumps(report_json, indent=4)
            
            f.write("var report_data = " + json_str.replace('Nan', '""').replace('NaN', '""'))
    
    def png_dict(self, png_dic):
        png_json = dict()
        id2key = [
            ("102002", "structure_snp_stat"),
            ("1018", "structure_rmat"),
            ("1019", "structure_fusion"),
            ("1014", "geneset_enrich_do"),
            ("1015", "allgene_venn"),
            ("1016", "allgene_heat"),
            ("1017", "allgene_pca"),
            ("1010", "geneset_class_do"),
            ("1011", "geneset_enrich_go"),
            ("1012", "geneset_enrich_kegg"),
            ("1013", "geneset_enrich_reactome"),
            ("101802", "structure_rmat_diffs"),
            ("1025", "alignments_region"),
            ("1024", "alignments_coverage"),
            ("1026", "allgene_distribution"),
            ("1021", "qc_qual_err"),
            ("1020", "structure_snp"),
            ("1023", "qc_qual_box"),
            ("1022", "qc_base_line"),
            ("1009", "geneset_class_reactome"),
            ("1008", "geneset_class_kegg"),
            ("1007", "geneset_class_go"),
            ("1005", "geneset_heatmap"),
            ("1004", "diff_venn"),
            ("1003", "diff_ma"),
            ("1002", "diff_volcano"),
            ("1001", "diff_summary")
        ]

        for id, k in id2key:
            if id in png_dic:
                png_json[k] = png_dic[id]["img"].replace(".png", ".jpg")
        return png_json
            

    def chart_json_batch(self, chart_json):
        report_json = {}
        overview = {
            "rawdata_total":"",
            "cleandata_min": "",
            "cleandata_q30_min":"",
            "mapping_min": "",
            "exp_num": "",
            "exp_known": "",
            "exp_novel": "",
            "exp_t_num": "",
            "exp_t_known": "",
            "exp_t_novel": "",
            
        }
        if "group_dict" in chart_json:
            table_header = [
                {"field": "group", "label": "group"},
                {"field": "sample", "label": "sample"}
            ]
            table = []
            for g,ss in chart_json["group_dict"].items():
                for s in ss:
                    table.append({
                        "sample": s,
                        "group": g
                    })
            report_json.update({"sample_data": {"table_header": table_header, "table": table}})

        if "raw_data_stat" in chart_json and "qc_data_stat" in chart_json:
            raw_df = pd.read_table(chart_json["raw_data_stat"], header=0, sep="\t")
            clean_df = pd.read_table(chart_json["qc_data_stat"], header=0, sep="\t")
            a = sum(raw_df["Total_Bases"]) / 1000000000.0
            a_str = str("%.2f" % a)
            b = min(clean_df["Total_Bases"]) / 1000000000.0
            b_str = str("%.2f" % b)
            c = min(clean_df["Q30%"])
            c_str = str(c) + "%"
            overview = {
                "rawdata_total":a_str,
                "cleandata_min": b_str,
                "cleandata_q30_min":c_str,
                "sample_num": len(raw_df)
            }



            table_header1, table1 = self.data_frame_to_json(raw_df)
            table_header2, table2 = self.data_frame_to_json(clean_df)
            raw_data = {
                "table_header": table_header1,
                "table": table1
            }
            clean_data = {
                "table_header": table_header2,
                "table": table2
            }

            report_json.update({"raw_data": raw_data, "clean_data": clean_data})

        if "align_stat" in chart_json:
            align_df = pd.read_table(chart_json["align_stat"], header=0, sep="\t")
            table_header, table = self.data_frame_to_json(align_df)
            align_data = {
                "table_header": table_header,
                "table": table
            }

            map_rate =  align_df["Total_mapped"].map(lambda x:x.split("(")[1].strip(")"))
            overview.update({
                "mapping_min": min(map_rate),
                "mapping_max": max(map_rate),
                
            })

            report_json.update({"align_data": align_data})

        if "gene_exp_all" in chart_json:
            g_a = pd.read_table(chart_json["gene_exp_all"], header=0, sep="\t", index_col=0)
            g_a_sum = g_a.apply(sum, axis="columns")
            g_a_num = len(g_a_sum[g_a_sum!=0])
        else:
            g_a_num = 0

        if "gene_exp_ref" in chart_json:
            g_r = pd.read_table(chart_json["gene_exp_ref"], header=0, sep="\t", index_col=0)
            g_r_sum = g_r.apply(sum, axis="columns")
            g_r_num = len(g_r_sum[g_r_sum!=0])
        else:
            g_r_num = 0
        if "tran_exp_all" in chart_json:
            t_a = pd.read_table(chart_json["tran_exp_all"], header=0, sep="\t", index_col=0)
            t_a_sum = t_a.apply(sum, axis="columns")
            t_a_num = len(t_a_sum[t_a_sum!=0])
        else:
            t_a_num = 0

        if "tran_exp_ref" in chart_json:
            t_r = pd.read_table(chart_json["tran_exp_ref"], header=0, sep="\t", index_col=0)
            t_r_sum = t_r.apply(sum, axis="columns")
            t_r_num = len(t_r_sum[t_r_sum!=0])
        else:
            t_r_num = 0

        overview.update({
            "exp_num": g_a_num,
            "exp_known": g_r_num,
            "exp_novel": g_a_num - g_r_num,
            "exp_t_num": t_a_num,
            "exp_t_known": t_r_num,
            "exp_t_novel": t_a_num -t_r_num
        })

        report_json.update({"overview": overview})


        if "subcluster_list" in chart_json:
            table_header = [
                {"field": "cluster", "label": "Cluster id"}, 
                {"field": "members", "Label": "Members"}
            ]
            table = []
            for f in chart_json["subcluster_list"]:
                [cluster, num] = f.split("subcluster_")[1].split(".")[0].split("_")
                table.append({"cluster":cluster, "members":num})

            geneset_cluster_data = {
                "table_header": table_header,
                "table": table
            }
            report_json.update({"geneset_cluster_data": geneset_cluster_data})

        if "go_class" in chart_json:
            go_class_files = [chart_json["go_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]] 
            go_class_df = pd.read_table(go_class_files[0], header=0, sep="\t", index_col=False)
            table_header, table = self.data_frame_to_json(go_class_df)
            geneset_go_class_data = {
                "table_header": table_header,
                "table": table[:10]
            }
            report_json.update({"geneset_class_go_data": geneset_go_class_data})

        if "do_class" in chart_json:
            do_class_files = [chart_json["do_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]] 
            do_class_df = pd.read_table(do_class_files[0], header=0, sep="\t", index_col=False)
            table_header, table = self.data_frame_to_json(do_class_df)
            geneset_do_class_data = {
                "table_header": table_header,
                "table": table[:10]
            }
            report_json.update({"geneset_class_do_data": geneset_do_class_data})

        if "kegg_class" in chart_json:
            kegg_class_files = [chart_json["kegg_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]] 
            kegg_class_df = pd.read_table(kegg_class_files[0], header=0, sep="\t", index_col=False)
            kegg_class_df["Ko_ids"] = kegg_class_df["Ko_ids"].map(lambda x: x[:50].replace("|", " ") + "...")
            
            table_header, table = self.data_frame_to_json(kegg_class_df)
            geneset_kegg_class_data = {
                "table_header": table_header,
                "table": table[:10]
            }
            report_json.update({"geneset_class_kegg_data": geneset_kegg_class_data})

        if "reactome_class" in chart_json:
            reactome_class_files = [chart_json["reactome_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]] 
            reactome_class_df = pd.read_table(reactome_class_files[0], header=0, sep="\t", index_col=False)
            table_header, table = self.data_frame_to_json(reactome_class_df)
            geneset_reactome_class_data = {
                "table_header": table_header,
                "table": table[:10]
            }

            report_json.update({"geneset_class_reactome_data": geneset_reactome_class_data})

        if "go_enrich" in chart_json:
            go_enrich_files = [chart_json["go_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]] 
            go_enrich_df = pd.read_table(go_enrich_files[0], header=0, sep="\t", index_col=False)
            go_enrich_df["seq_list"] = go_enrich_df["seq_list"].map(lambda x: x[:50].replace("|", " ") + "...")
            table_header, table = self.data_frame_to_json(go_enrich_df)
            geneset_go_enrich_data = {
                "table_header": table_header,
                "table": table[:10]
            }
            report_json.update({"geneset_enrich_go_data": geneset_go_enrich_data})

        if "do_enrich" in chart_json:
            do_enrich_files = [chart_json["do_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]] 
            do_enrich_df = pd.read_table(do_enrich_files[0], header=0, sep="\t", index_col=False)
            do_enrich_df["Genes"] = do_enrich_df["Genes"].map(lambda x: x[:50].replace("|", " ") + "...")
            # do_enrich_df["Genes_name"] = do_enrich_df["Genes_name"].map(lambda x: x[:50].replace("|", " ") + "...")

            table_header, table = self.data_frame_to_json(do_enrich_df)
            geneset_do_enrich_data = {
                "table_header": table_header,
                "table": table[:10]
            }
            report_json.update({"geneset_enrich_do_data": geneset_do_enrich_data})

        if "kegg_enrich" in chart_json:
            kegg_enrich_files = [chart_json["kegg_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]] 
            kegg_enrich_df = pd.read_table(kegg_enrich_files[0], header=0, sep="\t", index_col=False)
            kegg_enrich_df["Genes"] = kegg_enrich_df["Genes"].map(lambda x: x[:50].replace("|", " ") + "...")
            table_header, table = self.data_frame_to_json(kegg_enrich_df)
            geneset_kegg_enrich_data = {
                "table_header": table_header,
                "table": table[:10]
            }
            report_json.update({"geneset_enrich_kegg_data": geneset_kegg_enrich_data})

        if "reactome_enrich" in chart_json:
            reactome_enrich_files = [chart_json["reactome_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]] 
            reactome_enrich_df = pd.read_table(reactome_enrich_files[0], header=0, sep="\t", index_col=False)
            reactome_enrich_df["Genes"] = reactome_enrich_df["Genes"].map(lambda x: x[:50].replace("|", " ") + "...")
            table_header, table = self.data_frame_to_json(reactome_enrich_df)
            geneset_reactome_enrich_data = {
                "table_header": table_header,
                "table": table[:10]
            }

            report_json.update({"geneset_enrich_reactome_data": geneset_reactome_enrich_data})

        if "annot_detail" in chart_json:
            data_df = pd.read_table(chart_json["annot_detail"], header=0, sep="\t", index_col=False)
            clean_df = data_df[:10]

            table_header, table = self.data_frame_to_json(clean_df)
            annot_detail = {
                "table_header": table_header,
                "table": table
            }
            report_json.update({"annot_detail": annot_detail})


        if "diff_exp_summary" in chart_json:
            summary_pd = pd.read_table(chart_json["diff_exp_summary"], header=0, sep="\t", index_col=False)

            categories = sorted(list(summary_pd.columns[1:-1]))
            up_list = list()
            down_list = list()

            table_header = [
                {"field": "diff", "label": "diff_group"}, 
                {"field": "total", "label": "total DEG"},
                {"field": "up", "label": "up"},
                {"field": "down", "label": "down"},
            ]
            table = []
            for cmp1 in sorted(categories):
                stat = dict(summary_pd[cmp1].value_counts())
                up = stat.get("yes|up", 0)
                down = stat.get("yes|down", 0)
                table.append({
                    "diff": cmp1,
                    "total": up + down,
                    "up": up,
                    "down": down
                })

            diff_summary = {
                "table_header": table_header,
                "table": table
            }
            report_json.update({"diff_summary": diff_summary})

        return report_json



if __name__ == '__main__':
    a = DataReport()
    png_json = sys.argv[2]
    chart_json = sys.argv[1]
    a.run(chart_json, png_json)