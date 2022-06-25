# coding=utf-8
import xml.etree.ElementTree as ET
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
import shutil
import numpy as np
from scipy import stats
import xml.etree.ElementTree as ET
import lxml.html
from collections import OrderedDict
from collections import defaultdict

reload(sys)
sys.setdefaultencoding('utf-8')

class Chart(object):
    def __init__(self):
        """
        设置数据库，连接到mongod数据库，kegg_ko,kegg_gene,kegg_pathway_png三个collections2
        """

        """
        class: sample[A1, B1], group[A, B], classify: [samples, groups], level[G, T], kind[ref, new, all], category[mRNA, lncRNA]
        """
        self.sample_list = list()
        self.command_list = list()
        # chart_dir = "/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/work/sg_chart/test2"
        chart_dir = Config().SOFTWARE_DIR + "/bioinfo/sg_chart"
        self.phantomjs_dir = Config().SOFTWARE_DIR + "/program/phantomjs/phantomjs-2.1.1-linux-x86_64/bin/"
        self.wkhtml_dir = Config().SOFTWARE_DIR + "/install_packages/wkhtmltox-0.12.6-1/usr/local/bin/"
        self.mode_dir = chart_dir + "/refrna_v2"
        self.mode_mode = chart_dir + "/model"
        self.mode_mode2 = chart_dir + "/highchart_model"
        self.js_list = list()
        self.work_dir = os.environ.get("PWD", "") + "/"
        self.cairo_svg = Config().SOFTWARE_DIR + "/bioinfo/rna/miniconda3/bin/cairosvg"
        try:
            os.link(self.mode_mode + '/iconfont.woff', self.work_dir + '/iconfont.woff')
        except:
            pass

    def chart_raw_qc(self, sample_list, raw_qc_file_list, qc_type="raw"):
        for sample, qc_file in zip(sample_list, raw_qc_file_list):
            atcgn_pct = list()
            qc_file_list = qc_file.split(",")
            if len(qc_file_list) == 2:
                r_qc = pd.read_table(qc_file_list[0], header=0)
                r_qc["A_pct"] = r_qc["A_Count"]/r_qc["Max_count"] *100
                r_qc["T_pct"] = r_qc["T_Count"]/r_qc["Max_count"] *100
                r_qc["G_pct"] = r_qc["G_Count"]/r_qc["Max_count"] *100
                r_qc["C_pct"] = r_qc["C_Count"]/r_qc["Max_count"] *100
                r_qc["N_pct"] = r_qc["N_Count"]/r_qc["Max_count"] *100
                r_qc["err"] =  10**(r_qc["mean"]/-10) *100

                l_qc = pd.read_table(qc_file_list[1], header=0)
                l_qc["A_pct"] = l_qc["A_Count"]/l_qc["Max_count"] *100
                l_qc["T_pct"] = l_qc["T_Count"]/l_qc["Max_count"] *100
                l_qc["G_pct"] = l_qc["G_Count"]/l_qc["Max_count"] *100
                l_qc["C_pct"] = l_qc["C_Count"]/l_qc["Max_count"] *100
                l_qc["N_pct"] = l_qc["N_Count"]/l_qc["Max_count"] *100
                l_qc["err"] =  10**(l_qc["mean"]/-10) *100

                source_base = [
                    {
                        "data": list(r_qc[x + "_pct"]) + list(l_qc[x + "_pct"]),
                        "name": x
                    }
                    for x in ["A", "T", "C", "G", "N"]
                ]

                source_err = [
                    {
                        "data": list(r_qc["err"]) + list(l_qc["err"]),
                        "name": sample
                    }
                ]
                source_box = list()
                pos = 1
                for line_dict in l_qc.iterrows():
                    source_box.append({
                            "category": "All",
                            "data": [
                                line_dict[1]["lW"],
                                line_dict[1]["Q1"],
                                line_dict[1]["med"],
                                line_dict[1]["Q3"],
                                line_dict[1]["rW"],
                            ],
                            "name": pos
                        })
                    pos += 1

                for line_dict in r_qc.iterrows():
                    source_box.append({
                            "category": "All",
                            "data": [
                                line_dict[1]["lW"],
                                line_dict[1]["Q1"],
                                line_dict[1]["med"],
                                line_dict[1]["Q3"],
                                line_dict[1]["rW"],
                            ],
                            "name": pos
                        })
                    pos += 1
                categories = range(1, len(r_qc) + len(l_qc) + 1)
            else:
                r_qc = pd.read_table(qc_file_list[0], header=0)
                r_qc["A_pct"] = r_qc["A_Count"]/r_qc["Max_count"] *100
                r_qc["T_pct"] = r_qc["T_Count"]/r_qc["Max_count"] *100
                r_qc["G_pct"] = r_qc["G_Count"]/r_qc["Max_count"] *100
                r_qc["C_pct"] = r_qc["C_Count"]/r_qc["Max_count"] *100
                r_qc["N_pct"] = r_qc["N_Count"]/r_qc["Max_count"] *100
                r_qc["err"] =  10**(r_qc["mean"]/-10) *100
                source_base = [
                    {
                        "data": list(r_qc[x + "_pct"]),
                        "name": x
                    }
                    for x in ["A", "T", "C", "G", "N"]
                ]
                source_err = [
                    {
                        "data": list(r_qc["err"]),
                        "name": sample
                    }
                ]
                source_box = list()
                pos = 1
                for line_dict in r_qc.iterrows():
                    source_box.append({
                            "category": "All",
                            "data": [
                                line_dict[1]["lW"],
                                line_dict[1]["Q1"],
                                line_dict[1]["med"],
                                line_dict[1]["Q3"],
                                line_dict[1]["rW"],
                            ],
                            "name": pos
                        })
                    pos += 1
                categories = range(1, len(r_qc) + 1)

            self.chart_raw_qc_base(sample, ".{}_qc_base".format(qc_type), source_base, categories, "qc.base.line.json", qc_type=qc_type)
            self.chart_raw_qc_err(sample, ".{}_qc_error".format(qc_type), source_err, categories, "qc.err.line.json", qc_type=qc_type)
            self.chart_raw_qc_box(sample, ".{}_qc_qual".format(qc_type), source_box, "qc.qual.box.json", qc_type=qc_type)


    def chart_satu(self, sample_list, satu_file_list, curve_files):
        for sample, satu_file, curve_file in zip(sample_list, satu_file_list, curve_files):
            categaries = {"column1": "[0-0.3)", "column2": "[0.3-0.6)", "column3": "[0.6-3.5)", "column4": "[3.5-15)", "column5": "[15-60)", "column6": ">=60"}
            steps = map(str, range(5, 101, 5))
            with open(satu_file, "r") as f:
                for line in f:
                    if re.match(r"legend", line):
                        all_num = re.findall("num=[\d]*", line)
                        # print all_num
                        categaries = {
                            "column5": "[15-60)=" + all_num[4][4:],
                            "column4": "[3.5-15)=" + all_num[3][4:],
                            "column6": ">=60=" + all_num[5][4:],
                            "column1": "[0-0.3)=" + all_num[0][4:],
                            "column3": "[0.6-3.5)=" + all_num[2][4:],
                            "column2": "[0.3-0.6)=" + all_num[1][4:]
                        }
            source = list()
            source2 = [["name","x","y","category"]]
            with open(curve_file, "r") as f:
                n = 1
                for line in f:
                    pcts = line.strip("\n").split("\t")[1:]
                    source.append({
                        "data": map(float, pcts),
                        "group": categaries["column" + str(n)]
                    })
                    for step,pct in zip(steps, pcts):
                        source2.append([categaries["column" + str(n)], int(step), float(pct), categaries["column" + str(n)]])
                    n += 1

            steps = map(str, range(5, 101, 5))
            self.chart_satu_line(sample, ".align.satu", source, source2, steps, "align.satu.line.json")


    def chart_raw_qc_base(self, name, out, source, categories, json_mode, qc_type=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            if qc_type:
                a["title"]["text"] = a["title"]["text"].replace("clean", qc_type)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([name + out + ".line.js", {}])

    def chart_raw_qc_err(self, name, out, source, categories, json_mode, qc_type=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            if qc_type:
                a["title"]["text"] = a["title"]["text"].replace("clean", qc_type)
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {}])

    def chart_raw_qc_box(self, name, out, source, json_mode, qc_type=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".box.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            if qc_type:
                a["title"]["text"] = a["title"]["text"].replace("clean", qc_type)
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".box.js", {}])

    def chart_satu_line(self, name, out, source, source2, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            a["dataset"][1]["source"] = source2
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {}])

    def chart_line_point(self, name, out, source, source_point, categories, title, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = title
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            a["dataset"][1]["source"] = source_point
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {}])


    def chart_coverage(self, sample_list, coverage_file_list):
        steps = map(str, range(1, 101, 1))
        source_all = list()
        for sample, cov_file in zip(sample_list, coverage_file_list):
            source = list()
            with open(cov_file, "r") as f:
                f.readline()
                line = f.readline()
                depths = line.strip("\n").split("\t")[1:]
                source.append({
                    "data": map(int, map(float, depths)),
                    "name": sample
                })
                source_all.append({
                    "data": map(int, map(float, depths)),
                    "name": sample
                })
                self.chart_coverage_line(sample, ".align_coverage", source, steps, "align.coverage.line.json")

        self.chart_coverage_line("all", ".align_coverage", source_all, steps, "align.coverage.line.json")

    def chart_coverage_line(self, name, out, source, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {}])


    def chart_readpos(self, sample_list, readpos_file_list):
        for sample, readpos_file in zip(sample_list, readpos_file_list):
            source = [["value","name","category"]]
            values = list()
            with open(readpos_file, "r") as f:
                for line in f.readlines()[5:-1]:
                    line = line.strip().split()
                    values.append(float(line[2]))
            values_new = values[:4]
            values_new.append(sum([values[6],values[9]]))
            total = sum(values_new)
            source += [
                [values_new[3], "Introns", "Introns:{}%".format("%0.4f" %(values_new[3]/total * 100))],
                [values_new[2], "3'UTR", "3'UTR:{}%".format("%0.4f" %(values_new[2]/total * 100))],
                [values_new[0], "CDS", "CDS:{}%".format("%0.4f" %(values_new[0]/total * 100))],
                [values_new[1], "5'UTR","5'UTR:{}%".format("%0.4f" %(values_new[1]/total * 100))],
                [values_new[4], "Intergenic", "Intergenic:{}%".format("%0.4f" %(values_new[4]/total * 100))]
            ]

            title = "Percent  of reads mapped to genome regions ({sample_name})".format(sample_name = sample)
            self.chart_pie(sample, ".align_pos_dist", source, "align.pos_dist.pie.json", title)

    def chart_pie(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".pie.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".pie.js", {}])

    def chart_readchr(self, sample_list, readchr_file_list):
        for sample, readchr_file in zip(sample_list, readchr_file_list):
            read_chr = pd.read_table(readchr_file, header=0)
            read_chr_sorted = read_chr.sort_values(by="read_num", ascending=False)[:30]
            source = [{
                "data": list(read_chr_sorted["read_num"]),
                "name": sample
            }]
            categories = list(read_chr_sorted["#chr"])

            self.chart_readchr_line(sample, ".align_chr_dist", source, categories, "align.chr_dist.bar.json")

    def chart_readchr_line(self, name, out, source, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".bar.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories

            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".bar.js", {}])

    def chart_bar_and_line(self, name, out, source_bar, source_line, categories_line, geneset_name, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".bar_line.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(geneset_name = geneset_name)
            a["dataset"][0]["source"] = source_bar
            a["dataset"][1]["source"] = source_line
            a["dataset"][1]["categories"] = categories_line
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".bar_line.js", {}])

    def chart_circ(self, name, out, source_circ, source_class, title, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".circ.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = title
            a["dataset"][0]["source"] = source_circ
            a["dataset"][1]["source"] = source_class
            a["series"][0]["visualMap"][0]["data"] = len(source_class) + 1
            a['legend'][1]['title']['text'] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".circ.js", {}])

    def chart_dag(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".dag.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            a["dataset"][1]["source"] = source
            a["chart"]["width"] = 1800
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".dag.js", {}])

    def chart_network(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".network.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"][0] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".network.js", {"type": "network"}])

    def chart_bar(self, name, out, source_bar, geneset_name, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".bar.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(geneset_name = geneset_name)
            a["dataset"][0]["source"] = source_bar
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".bar.js", {"delay": 2000}])

    def chart_buble(self, name, out, source_buble, geneset_name, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".buble.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(geneset_name = geneset_name)
            a["dataset"][0]["source"] = source_buble
            size_list = [s[2] for s in source_buble[1:]]
            if len(size_list) > 0:
                a["legend"][1]["customizeSize"] = [int(min(size_list)),
                                                   int(min(size_list) + 0.33 * (max(size_list) - min(size_list)) + 0.5),
                                                   int(min(size_list) + 0.66 * (max(size_list) - min(size_list)) + 0.5),
                                                   int(max(size_list))]
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".buble.js", {"width": 900}])

    def chart_buble2(self, name, out, source_buble_list, geneset_name, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".buble2.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(geneset_name = geneset_name)
            a["dataset"] = [{"rawName": True, "source": source_buble, "colName": True, "sourceType": "matrix"}
                            for source_buble in source_buble_list]
            size_list = []
            for source_buble in source_buble_list:
                for s in source_buble[1:]:
                    size_list.append(s[2])
            if len(size_list) > 0:
                a["legend"][1]["customizeSize"] = [int(min(size_list)),
                                                   int(min(size_list) + 0.33 * (max(size_list) - min(size_list)) + 0.5),
                                                   int(min(size_list) + 0.66 * (max(size_list) - min(size_list)) + 0.5),
                                                   int(max(size_list))]

            # 分类不足时减少坐标轴
            a["yAxis"] = a["yAxis"][:len(source_buble_list)]
            a["xAxis"] = a["xAxis"][:len(source_buble_list)]
            a["series"] = a["series"][:len(source_buble_list)]

            for x in a["xAxis"]:
                x["length"] = "{}%".format(float(100)/len(source_buble_list))

            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".buble2.js", {}])

    def chart_assemble(self, assemble_step):
        # trans_count_stat_200.txt
        with open(assemble_step, "r") as f:
            f.readline()
            source = [["item"], ["series"], ["category"]]
            for line in f.readlines()[:-1]:
                cols = line.strip().split("\t")
                source[0].append(cols[0])
                source[1].append(int(cols[1]))
                source[2].append(cols[0])
        self.chart_column2("all", ".assemble_len", source, "assemble.len.column.json")

    def chart_assemble_relation(self, ref_gene2tran, ref_tran2exon, new_gene2tran, new_tran2exon):
        ref_g2t_pd = pd.read_table(ref_gene2tran, header=None)
        ref_g2t_source = [
            list(ref_g2t_pd[1])
        ]

        categories = list(ref_g2t_pd[0])
        title = "Transcripts per gene(ref)"
        self.chart_columns_highchart("ref", ".assemble_relation_g2t", ref_g2t_source, categories, "assemble.gene_trans.column.json", title)

        ref_t2e_pd = pd.read_table(ref_gene2tran, header=None)
        ref_t2e_source = [
            list(ref_t2e_pd[1])
        ]
        categories = list(ref_t2e_pd[0])
        title = "Exons per transcript(ref)"
        self.chart_line_highchart("ref", ".assemble_relation_t2e", ref_t2e_source, categories, "assemble.trans_exon.line.json", title)

        if new_gene2tran:
            new_g2t_pd = pd.read_table(new_gene2tran, header=None)
            new_g2t_source = [
                list(new_g2t_pd[1])
            ]

            categories = list(new_g2t_pd[0])
            title = "Transcripts per gene(new)"
            self.chart_columns_highchart("new", ".assemble_relation_g2t", new_g2t_source, categories, "assemble.gene_trans.column.json", title)

            new_t2e_pd = pd.read_table(new_gene2tran, header=None)
            new_t2e_source = [
                list(new_t2e_pd[1])
            ]
            categories = list(new_t2e_pd[0])
            title = "Exons per transcript(new)"
            self.chart_line_highchart("new", ".assemble_relation_t2e", new_t2e_source, categories, "assemble.trans_exon.line.json", title)


    def chart_column2(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".column.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".column.js", {}])

    def chart_column_and_conf(self, name, out, source, source2, color, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".column_conf.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            a["dataset"][1]["source"] = source2
            a["series"][0]["visualMap"][0]["visualColorValue"] = color
            a["legend"][0]["color"] = color
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".column_conf.js", {}])

    def chart_assemble_new(self, code_file):
        # trans_count_stat_200.txt
        with open(code_file, "r") as fr:
            source =  [["name","value","category"]]
            code_list = list()
            num_list = list()
            for line in fr:
                new_gene_list = []
                lines = line.strip().split("\t")
                gene_list = lines[1].strip().split(",")
                code_list.append(lines[0])
                num_list.append(len(gene_list))

            total = sum(num_list)
            for c,n in zip(code_list, num_list):
                source.append([
                    c,
                    n,
                    "{}:{}%".format(c, "%0.4f" %(float(n)/total * 100))
                ])

        self.chart_assemble_new_pie("all", ".assemble_new", source, "assemble.new.pie.json")

    def chart_assemble_new_pie(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".pie.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".pie.js", {}])

    def chart_pie(self, name, out, source, json_mode, title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".pie.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            if title:
                a["title"]["text"] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".pie.js", {}])

    def get_exp_list(self, exp_file):
        exp_g_set = set()
        if os.path.exists(exp_file):
            with open(exp_file, 'r') as gene_exp_f:
                gene_exp_f.readline()
                for line in gene_exp_f:
                    exps = line.strip().split("\t")[1:]
                    if sum(map(float, exps)) != 0:
                        exp_g_set.add(line.strip().split("\t")[0])
        return exp_g_set


    def get_venn_list(self, venn_dir, db, level="G"):
        tail = "gene"

        database_venn = {
            'NR': 'nr/nr_venn_{}.txt'.format(tail),
            'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
            'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
            'GO': 'go/go_venn_{}.txt'.format(tail),
            'COG': 'cog/cog_venn_{}.txt'.format(tail),
        }
        tail = "tran"
        database_venn_tran = {
            'NR': 'nr/nr_venn_{}.txt'.format(tail),
            'Swiss-Prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Swiss-prot': 'swissprot/swissprot_venn_{}.txt'.format(tail),
            'Pfam': 'pfam/pfam_venn_{}.txt'.format(tail),
            'KEGG': 'kegg/kegg_venn_{}.txt'.format(tail),
            'GO': 'go/go_venn_{}.txt'.format(tail),
            'COG': 'cog/cog_venn_{}.txt'.format(tail),
        }
        v3_db2db = {
            "pfam": 'Pfam',
            "kegg": "KEGG",
            "swissprot": "Swiss-Prot",
            "cog": "COG",
            "go": "GO",
            "nr": "NR",
            "annotation": "Total_anno",
            "total": "Total",
            "Pfam": 'Pfam',
            "KEGG": "KEGG",
            "Swiss-Prot": "Swiss-Prot",
            "COG": "COG",
            "GO": "GO",
            "NR": "NR",
            "Total_anno": "Total_anno",
            "Total": "Total"
        }
        if level == "G":
            venn_file = venn_dir + "/" + database_venn[db]
        else:
            venn_file = venn_dir + "/" + database_venn_tran[db]
        with open(venn_file, "rb") as f:
            venn_list = map(lambda x: x.strip("\n"), f.readlines())
        return venn_list


    def chart_annotation_stat(self, names, gene_exp, trans_exp, annot_stat, venn_dir):
        # trans_count_stat_200.txt
        exp_g_set = self.get_exp_list(gene_exp)
        exp_t_set = self.get_exp_list(trans_exp)

        # print exp_g_set

        source_gene_venn = []
        source_tran_venn = []

        with open(annot_stat, "r") as fr:
            source_gene =  [["item"], ["series"], ["category"]]
            source_tran =  [["item"], ["series"], ["category"]]
            source_gene_venn = []
            source_tran_venn = []

            source_gene_exp =  [["item"], ["series"], ["category"]]
            source_tran_exp =  [["item"], ["series"], ["category"]]
            source_gene_venn_exp = []
            source_tran_venn_exp = []

            for line in fr.readlines()[1:-2]:

                cols = line.strip().split("\t")
                # print cols
                source_tran[0].append(cols[0])
                source_tran[1].append(int(cols[1]))
                source_tran[2].append(cols[0])
                source_gene[0].append(cols[0])
                source_gene[1].append(int(cols[2]))
                source_gene[2].append(cols[0])
                gene_venn_list = self.get_venn_list(venn_dir, cols[0], level="G")
                tran_venn_list = self.get_venn_list(venn_dir, cols[0], level="T")
                source_gene_venn.append({
                    "data": gene_venn_list,
                    "name": cols[0]
                })
                source_tran_venn.append({
                    "data": tran_venn_list,
                    "name": cols[0]
                })

                gene_venn_list_exp = set(exp_g_set) & set(gene_venn_list)
                tran_venn_list_exp = set(exp_t_set) & set(tran_venn_list)
                # print gene_venn_list_exp
                source_gene_venn_exp.append({
                    "data": list(gene_venn_list_exp),
                    "name": cols[0]
                })
                source_tran_venn_exp.append({
                    "data": list(tran_venn_list_exp),
                    "name": cols[0]
                })

                source_tran_exp[0].append(cols[0])
                source_tran_exp[1].append(len(tran_venn_list_exp))
                source_tran_exp[2].append(cols[0])
                source_gene_exp[0].append(cols[0])
                source_gene_exp[1].append(len(tran_venn_list_exp))
                source_gene_exp[2].append(cols[0])


        self.chart_annotation_stat_column("", "annot_gene_stat", source_gene, "annotation.stat.column.json")
        self.chart_annotation_stat_column("", "annot_tran_stat", source_tran, "annotation.stat.column.json")
        self.chart_annotation_stat_column("", "annot_gene_stat_exp", source_gene_exp, "annotation.stat.column.json")
        self.chart_annotation_stat_column("", "annot_tran_stat_exp", source_tran_exp, "annotation.stat.column.json")
        self.chart_annotation_stat_venn("", "annot_gene_stat_venn", source_gene_venn, "annotation.stat.venn.json")
        self.chart_annotation_stat_venn("", "annot_tran_stat_venn", source_tran_venn, "annotation.stat.venn.json")
        self.chart_annotation_stat_venn("", "annot_gene_stat_venn_exp", source_gene_venn_exp, "annotation.stat.venn.json")
        self.chart_annotation_stat_venn("", "annot_tran_stat_venn_exp", source_tran_venn_exp, "annotation.stat.venn.json")
        self.chart_annotation_stat_upset("", "annot_gene_stat_upset", source_gene_venn, "annotation.stat.upset.json")
        self.chart_annotation_stat_upset("", "annot_tran_stat_upset", source_tran_venn, "annotation.stat.upset.json")
        self.chart_annotation_stat_upset("", "annot_gene_stat_upset_exp", source_gene_venn_exp, "annotation.stat.upset.json")
        self.chart_annotation_stat_upset("", "annot_tran_stat_upset_exp", source_tran_venn_exp, "annotation.stat.upset.json")


    def chart_annotation_stat_column(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".column.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".column.js", {}])

    def chart_annotation_stat_venn(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".venn.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".venn.js", {"convert": "wkhtmltopdf", "zoom": 7, "width": 720, "height": 600}])

    def chart_annotation_stat_upset(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".upset.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".upset.js", {}])

    def chart_venn(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".venn.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".venn.js", {"convert": "wkhtmltopdf", "zoom": 7, "width": 900, "height": 600}])

    def chart_upset(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".upset.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".upset.js", {}])


    def get_box(self, all_exp_pd):
        """
        get box plot info for each column of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        """
        stat_dict_list = list()
        target_columns = all_exp_pd.columns
        for each in target_columns:
            exp_pd = all_exp_pd[each]
            exp_pd = exp_pd[exp_pd != 0]
            summary = exp_pd.describe()
            summary.index = [u'count', u'mean', u'std', u'min', u'q1', u'median', u'q3', u'max']
            tmp_dict = summary.to_dict()
            lt25 = exp_pd[exp_pd <= tmp_dict['q1']].shape[0]
            lt50 = exp_pd[exp_pd <= tmp_dict['median']].shape[0]
            lt75 = exp_pd[exp_pd <= tmp_dict['q3']].shape[0]
            upper_whisker = tmp_dict['q3'] + 1.5 * (tmp_dict['q3'] - tmp_dict['q1'])
            lower_whisker = tmp_dict['q1'] - 1.5 * (tmp_dict['q3'] - tmp_dict['q1'])
            upper_outliers = list(exp_pd[exp_pd > upper_whisker])
            lower_outliers = list(exp_pd[exp_pd < lower_whisker])
            tmp_dict.update({
                'sample': each,
                'upper_whisker': upper_whisker,
                'lower_whisker': lower_whisker,
                'upper_outliers': upper_outliers,
                'lower_outliers': lower_outliers,
            })
            stat_dict_list.append(tmp_dict)
        return stat_dict_list

    def process_exp_matrix(self, exp_matrix, log_base=None, group_dict=None):
        if type(exp_matrix) == str or type(exp_matrix) == bytes or isinstance(exp_matrix, unicode):
            all_exp_pd = pd.read_table(exp_matrix, index_col=0, header=0)
        else:
            print(exp_matrix, 'is assumed to be a pandas DataFrame Object')
            all_exp_pd = exp_matrix
        all_exp_pd.index.name = 'seq_id'

        if group_dict is not None:
            group_exp = list()
            for g in group_dict:
                g_exp = all_exp_pd.loc[:, group_dict[g]].mean(axis=1)
                g_exp.name = g
                group_exp.append(g_exp)
            all_exp_pd = pd.concat(group_exp, axis=1)
        all_exp_df = all_exp_pd.copy()
        if log_base:
            if log_base == math.e:
                all_exp_pd = np.log(all_exp_pd)
            elif log_base == 2:
                all_exp_pd = np.log2(all_exp_pd)
            elif log_base == 10:
                all_exp_pd = np.log10(all_exp_pd)
            else:
                print('log base of %s is not supported' %log_base)
        if len(all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]) < 2:
            all_exp_pd = np.log10(all_exp_df + 0.001)
        # return
        return all_exp_pd

    def get_box_source(self, exp_pd):
        # 获取box source
        stat_list = self.get_box(exp_pd)
        source_all = []
        for box in stat_list:
            source_all.append({
                "category": box['sample'],
                "data": [
                    box['lower_whisker'],
                    box['q1'],
                    box['median'],
                    box['q3'],
                    box['upper_whisker']
                ],
                "name": box['sample']
            })
        return source_all

    def get_density(self, all_exp_pd):
        """
        sampling 1000 density point for each columns of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        """
        records = list()
        target_columns = all_exp_pd.columns
        for sample in target_columns:
            exp = all_exp_pd[sample]
            exp = exp[exp != 0]
            density_func = stats.gaussian_kde(exp)
            min_exp, max_exp = exp.min(), exp.max()
            x_data = np.linspace(min_exp, max_exp, num=1000, endpoint=False)
            y_data = density_func(x_data)
            point_dict_list = {'log2exp': x_data, 'density': y_data}
            records.append(dict(sample=sample, data=point_dict_list))
        return records

    def get_density_source(self, exp_pd):
        # 获取box source
        stat_list = self.get_density(exp_pd)
        source_all = []
        for density in stat_list:
            source_all.append({
                "data": list(density["data"]["density"]), # 有点问题每个样本的横坐标不是一个体系，category如何设置
                "name": density['sample']
            })

        category = list(stat_list[0]["data"]["log2exp"])
        return category, source_all

    def get_violin_source(self, exp_pd):
        # 获取box source
        source_all = [["data", "name"]]
        target_columns = exp_pd.columns
        for sample in target_columns:
            exp = exp_pd[sample]
            source_all.append([list(exp), sample])
        return source_all

    def chart_exp_dis(self, names, gene_exp, tran_exp, group_dict):
        # 使用参考转录本做
        if os.path.exists(tran_exp):
            levels = ["G", "T"]
        else:
            levels = ["G"]
        for level in levels:
            for classify in ["samples", "groups"]:
                if classify == "samples":
                    g = None
                else:
                    g = group_dict
                if level == "G":
                    exp = gene_exp
                else:
                    exp = tran_exp

                all_exp_pd = self.process_exp_matrix(exp, log_base=10, group_dict=g)

                all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
                source_all = self.get_box_source(all_exp_pd)
                self.chart_raw_exp_box("{}_{}".format(level, classify), ".exp_distribution", source_all, "exp.dist.box.json")

                category, source_density = self.get_density_source(all_exp_pd)
                self.chart_raw_exp_density("{}_{}".format(level, classify), ".exp_distribution", source_density, category, "exp.dist.density.json")

                all_exp_pd = all_exp_pd.sample(frac=0.7)
                # all_exp_pd.reset_index(level=0, inplace=True)
                source_violin = self.get_violin_source(all_exp_pd)
                self.chart_raw_exp_violin("{}_{}".format(level, classify), ".exp_distribution", source_violin, "exp.dist.violin.json", y_title=None)

    def chart_exp_dis_one(self, exp, group_dict, exp_type="TPM"):
        # 使用参考转录本做
        all_exp_pd = self.process_exp_matrix(exp, log_base=10, group_dict=group_dict)

        all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
        source_all = self.get_box_source(all_exp_pd)
        title = "log10({})".format(exp_type)
        self.chart_raw_exp_box("", "exp_distribution", source_all, "exp.dist.box.json", title)

        category, source_density = self.get_density_source(all_exp_pd)
        self.chart_raw_exp_density("", "exp_distribution", source_density, category, "exp.dist.density.json")

        all_exp_pd = all_exp_pd.sample(frac=0.7)
        # all_exp_pd.reset_index(level=0, inplace=True)
        source_violin = self.get_violin_source(all_exp_pd)
        self.chart_raw_exp_violin("", "exp_distribution", source_violin, "exp.dist.violin.json", title)

    def chart_raw_exp_box(self, name, out, source, json_mode, y_title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".box.js", 'w') as fo:
            a = json.loads(f.read())
            if y_title:
                a["yAxis"][0]["title"]["text"] = y_title
            a["dataset"][0]["source"] = source
            word_list = [source1["name"] for source1 in source]
            a=self.reset_margin(a, word_list = word_list)
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".box.js", {}])

    def chart_raw_exp_density(self, name, out, source, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".density.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".density.js", {}])

    def chart_raw_exp_violin(self, name, out, source, json_mode, y_title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".violin.js", 'w') as fo:
            a = json.loads(f.read())
            if y_title:
                a["yAxis"][0]["title"]["text"] = y_title
            a["dataset"][0]["source"] = source
            word_list = [source1[1] for source1 in source]
            a=self.reset_margin(a, word_list = word_list)
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".violin.js", {}])

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
        if len(source_venn) <=6 :
            self.chart_exp_upset("all", ".exp", source_venn, "annotation.stat.upset.json")

    def chart_exp_venn_venn(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".venn.js", 'w') as fo:
            a = json.loads(f.read())

            a["dataset"][0]["source"] = source
            geneset_list = [s['name'][10:] for s in source]
            a = self.reset_margin(j=a, margin_type="right", word_list=geneset_list)
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".venn.js", {"convert": "wkhtmltopdf", "zoom": 7, "width": 720, "height": 600}])

    def chart_exp_upset(self, name, out, source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".upset.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".upset.js", {}])

    def chart_exp_pca(self, exp_pca_file, exp_pca_var_file, group_dict=None, exp_pca_ellipse=None, pcs=["PC1", "PC2"]):
        sample2group = dict()
        if group_dict:
            for g,ss in group_dict.items():
                for s in ss:
                    sample2group[s] = g

        pca_source = [["name","x","y","value","category"]]
        pca_df = pd.read_table(exp_pca_file, header=0)
        for line_dict in pca_df.iterrows():
            pca_source.append([
                line_dict[1]['sample'],
                line_dict[1][pcs[0]],
                line_dict[1][pcs[1]],
                line_dict[1]['sample'],
                sample2group.get(line_dict[1]['sample'], "all")
            ])

        t = pd.read_table(exp_pca_var_file, header=None)
        pc_ratio_dict = OrderedDict(zip(t[0], t[1]))

        x_title = "{}({}%)".format(pcs[0], "%0.2f" %(float(pc_ratio_dict[pcs[0]]) * 100))
        y_title = "{}({}%)".format(pcs[0], "%0.2f" %(float(pc_ratio_dict[pcs[1]]) * 100))

        self.chart_exp_pca_scatter("all", ".exp_relation_pca", pca_source, x_title, y_title, "exp.relation_pca.scatter.json")

        if exp_pca_ellipse:
            ellipse_source = list()
            pca_df = pd.read_table(exp_pca_ellipse, header=0)
            pca_df = pca_df.loc[pca_df["PC"] == "_".join(pcs), ]
            for line_dict in pca_df.iterrows():
                e_dict = dict(line_dict[1][2:])
                e_dict["cx"] = e_dict["m1"]
                e_dict["cy"] = e_dict["m2"]
                e_dict.pop("m1")
                e_dict.pop("m2")
                ellipse_source.append({
                    "group": line_dict[1]["group"],
                    "name": line_dict[1]["group"],
                    "ellipse": dict(e_dict)
                })

            self.chart_exp_pca_scatter2("all", ".exp_relation_pca_ell", pca_source, ellipse_source, x_title, y_title, "exp.relation_pca.scatter2.json")

    def chart_busco(self, summary_dict, text_dict):
        persent_dict = {'Complete(C) and single-copy(S)': [],
                        'Complete(C) and duplicated(D)': [],
                        'Fragmented(F)': [],
                        'Missing(M)': []}
        sample_list = list()
        text_list = list()
        print summary_dict
        print text_dict
        for sample in summary_dict:
            sample_list.append(sample)
            text_list.append({'y': sample, 'text': text_dict[sample]})
            summary = summary_dict[sample]
            with open(summary, 'r') as s:
                lines = s.readlines()[2:-1]
                persent_dict['Complete(C) and single-copy(S)'].append(float(lines[0].strip().split('\t')[2]))
                persent_dict['Complete(C) and duplicated(D)'].append(float(lines[1].strip().split('\t')[2]))
                persent_dict['Fragmented(F)'].append(float(lines[2].strip().split('\t')[2]))
                persent_dict['Missing(M)'].append(float(lines[3].strip().split('\t')[2]))
        # with open(summary, 'r') as s:
        #     persent_list = list()
        #     lines = s.readlines()[2:-1]
        #     for line in lines:
        #         persent = line.strip().split('\t')[2]
        #         persent_list.append(persent)
        type_list = ['Complete(C) and single-copy(S)', 'Complete(C) and duplicated(D)', 'Fragmented(F)', 'Missing(M)']
        # a = zip(persent_list, type_list)
        source_0 = [{'y': persent_dict[i], 'type': i} for i in type_list]


        # source_0 = [{'y': [float(i[0])], 'type': i[1]} for i in a]
        category_0 = sample_list
        source_1 = text_list
        self.chart_busco_stackbar("busco", ".bar", source_0, category_0, source_1,
                                   "busco.json")

    def chart_busco_stackbar(self, name, out, source_0, category_0, source_1, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".busco.js", 'w') as fo:
            a = json.loads(f.read())
            a['dataset'][0]['source'] = source_0
            a['dataset'][0]['categories'] = category_0
            a['dataset'][1]['source'] = source_1
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".busco.js", {}])

    def chart_scatter(self, name, out, source, x_title, y_title, title=None, json_mode=None):
        # 旧版插件


        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".scatter.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            a["xAxis"][0]["text"] = x_title
            a["xAxis"][0]["title"]["text"] = x_title
            a["yAxis"][0]["text"] = y_title
            a["yAxis"][0]["title"]["text"] = y_title
            a["title"]["text"] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".scatter.js", {}])


    def chart_exp_pca_scatter(self, name, out, source, x_title, y_title, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".scatter.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source
            a["xAxis"][0]["text"] = x_title
            a["xAxis"][0]["title"]["text"] = x_title
            a["yAxis"][0]["text"] = y_title
            a["yAxis"][0]["title"]["text"] = y_title
            categories = set([s[-1] for s in source])
            if len(categories) > 5:
                symbol = ["circle", "triangle", "diamond", "square", "triangle-down"]
                a["series"][0]["visualMap"][0]["visualSymbolValue"] = [symbol[i%5] for i in range(0, len(categories))]
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".scatter.js", {}])

    def chart_exp_pca_scatter2(self, name, out, source, source_ellipse, x_title, y_title, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".scatter.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source_ellipse
            a["dataset"][1]["source"] = source
            a["xAxis"][0]["text"] = x_title
            a["xAxis"][0]["title"]["text"] = x_title
            a["yAxis"][0]["text"] = y_title
            a["yAxis"][0]["title"]["text"] = y_title
            categories = set([s[-1] for s in source])
            if len(categories) > 5:
                symbol = ["circle", "triangle", "diamond", "square", "triangle-down"]
                a["series"][0]["visualMap"][0]["visualSymbolValue"] = [symbol[i%5] for i in range(0, len(categories))]
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".scatter.js", {}])

    def chart_exp_corr(self, sample_tree, sample_corr_file, group_dict=None):
        with open(sample_tree, 'r') as f:
            sample_tree = f.readline().strip()
            samples = f.readline().strip().split(";")

        print samples
        corr_result = pd.read_table(sample_corr_file, index_col=0, header=0)
        corr_result = corr_result.loc[samples, :]
        corr_source = [[""] + samples]
        for line_dict in corr_result.iterrows():
            corr_source.append(
                [line_dict[0]] + [line_dict[1][s] for s in samples]
            )
        tree_source = sample_tree

        sample2group_source = [["name", "group"]]
        if group_dict:
            for g,ss in group_dict.items():
                for s in ss:
                    sample2group_source.append([s, g])
        else:
            for s in samples:
                sample2group_source.append([s, s])

        self.chart_heat_tree("all", ".exp", corr_source , sample_tree, sample_tree, sample2group_source, sample2group_source, "exp.relation.heat_tree.json")

    def chart_heat_tree(self, name, out, corr_source, sample_tree, gene_tree, sample2group_source, gene2group_source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".heat_corr.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = corr_source
            a["dataset"][1]["source"] = gene2group_source
            if gene_tree:
                a["dataset"][1]["categories"] = gene_tree
            else:
                del a["dataset"][1]["categories"]
            a["dataset"][2]["source"] = sample2group_source
            if sample_tree:
                a["dataset"][2]["categories"] = sample_tree
            else:
                del a["dataset"][2]["categories"]
                a["series"] = a["series"][:2]
                a["dataset"] = a["dataset"][:2]
                a["legend"] = a["legend"][:2]

            if len(gene2group_source) <= 1:
                ##无基因聚类
                if sample_tree:
                    a["dataset"] = [a["dataset"][0], a["dataset"][2]]
                    a["legend"] = [a["legend"][0],  a["legend"][2]]
                    a["series"] = [a["series"][0], a['series'][2]]
                    a["series"][1]["datasetIndex"] = 1
                    a["series"][1]["visualMap"][0]["legendIndex"] = 1
                    a["xAxis"][0]["length"] = "95%"
                else:
                    a["dataset"] = [a["dataset"][0]]
                    a["legend"] = [a["legend"][0]]
                    a["series"] = [a["series"][0]]
                    a["xAxis"][0]["length"] = "95%"

            # 类型待确定
            a["chart"]["type"] = "heatmap_tree"
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".heat_corr.js", {"to_type": "svg"}])

    def chart_heat(self, name, out, corr_source, group_source, p_source, color, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".heat_corr.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = corr_source
            a["dataset"][1]["source"] = group_source
            a["dataset"][2]["source"] = p_source
            a["series"][1]["visualMap"][0]["visualValue"] = color
            a["series"][1]["visualMap"][0]["visualColorValue"] = color
            a["legend"][1]["color"] = color
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".heat_corr.js", {}])

    def chart_tree(self, name, out, sample_tree, sample2group_source, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".heat_corr.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = sample2group_source
            a["dataset"][0]["categories"] = sample_tree

            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".heat_corr.js", {}])


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
            else:
                up_list.append(0)
            if "yes|down" in stat:
                down_list.append(stat["yes|down"])
            else:
                down_list.append(0)
        all_source = [
            {
                "data": up_list,
                "name": "up"
            },
            {
                "data": down_list,
                "name": "down"
            }
        ]
        # 堆叠图
        all_source2 = [
            {
                "value": up_list,
                "category": "up"
            },
            {
                "value": down_list,
                "category": "down"
            }
        ]
        self.chart_diffexp_stat_bar("all", ".diffexp_summary", all_source, categories, "diffexp_stat.bar1.json")
        self.chart_diffexp_stat_bar2("all", ".diffexp_summary", all_source2, categories, "diffexp_stat.bar2.json")


    def reset_margin(self, j=None, margin_type="bottom", word_list=None):
        len_list = [len(x) for x in word_list]
        max_len = max(len_list)
        if max_len > 5:
            if "margin" in j["chart"]:
                if margin_type in j["chart"]["margin"]:
                    j["chart"]["margin"][margin_type] += (max_len - 5)*5
            else:
                j["chart"]["margin"] = {
                    margin_type: 30 + (max_len - 5)*5
                }
        return j


    def chart_diffexp_stat_bar(self, name, out, source, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".bar.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories

            a = self.reset_margin(j=a, margin_type="bottom", word_list=categories)
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".bar.js", {"delay": 2000}])

    def chart_diffexp_stat_bar2(self, name, out, source, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".bar2.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            a = self.reset_margin(j=a, margin_type="bottom", word_list=categories)
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".bar2.js", {}])

    def _get_volcano_status_cutoff(self, diff_table, pvalue_padjust):
        sig_status = list()
        sig_mark = diff_table['significant']
        reg_list = diff_table['regulate']
        if 'no' in list(sig_mark):
            no_sig_num = sig_mark[sig_mark == "no"].shape[0]
            sig_status.append('nosig_' + str(no_sig_num))
        if 'yes' in list(sig_mark):
            reg_mark = reg_list[sig_mark == 'yes']
            if 'down' in list(reg_mark):
                down_num = reg_mark[reg_mark == 'down'].shape[0]
                sig_status.append('down_' + str(down_num))
            if 'up' in list(reg_mark):
                up_num = reg_mark[reg_mark == 'up'].shape[0]
                sig_status.append('up_' + str(up_num))

        sig_pvalues = diff_table[pvalue_padjust][diff_table['significant'] == "yes"]
        log10_sig_pvalues = -np.log10(sig_pvalues)
        log10_pvalue_list = sorted(list(log10_sig_pvalues[log10_sig_pvalues > 0]))

        if len(sig_pvalues) > 2000:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.85)]
        elif len(sig_pvalues) > 1000:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.90)]
        elif len(sig_pvalues) > 500:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.95)]
        elif len(sig_pvalues) > 250:
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.99)]
        elif len(sig_pvalues) == 0:
            tmp = -np.log10(diff_table[pvalue_padjust])
            tmp_list = sorted(tmp[tmp > 0])
            if len(tmp_list) == 0:
                log10_pvalue_cutoff = 200
            else:
                log10_pvalue_cutoff = tmp_list[int(len(tmp_list) * 0.9)]
        else:
            # print(pvalue_padjust, diff_table, log10_pvalue_list)
            log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.8)]
        return sig_status, log10_pvalue_cutoff

    def _get_volcano_status_cutoff_noiseq(self, diff_table):
        sig_status = list()
        sig_mark = diff_table['significant']
        reg_list = diff_table['regulate']
        if 'no' in list(sig_mark):
            no_sig_num = sig_mark[sig_mark == "no"].shape[0]
            sig_status.append('nosig_' + str(no_sig_num))
        if 'yes' in list(sig_mark):
            reg_mark = reg_list[sig_mark == 'yes']
            if 'down' in list(reg_mark):
                down_num = reg_mark[reg_mark == 'down'].shape[0]
                sig_status.append('down_' + str(down_num))
            if 'up' in list(reg_mark):
                up_num = reg_mark[reg_mark == 'up'].shape[0]
                sig_status.append('up_' + str(up_num))

        return sig_status

    def chart_diffexp_scatter(self, diff_exp, cmps, soft=None, pvalue_padjust='padjust'):
        diff_pd = pd.read_table(diff_exp, header=0, sep='\t')
        if soft.lower() == "noiseq":
            need_cols = ['seq_id', 'fc', 'log2fc', 'D', 'prob', 'significant', 'regulate']
            volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', 'D', 'significant', 'regulate']]
            volcano_pd.columns = ['seq_id', 'log2fc', 'D', 'significant', 'regulate']
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
            volcano_pd_nosig.rename(columns={"D": "log10pvalue"}, inplace=True)
            volcano_pd_sig.rename(columns={"D": "log10pvalue"}, inplace=True)
        else:
            need_cols = ['seq_id', 'fc', 'log2fc', 'pvalue', 'padjust', 'significant', 'regulate']
            # pvalue_padjust='padjust'
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

        if soft.lower() == "noiseq":
            self.chart_diffexp_volcano(cmps, ".diffexp", data, title, "diffexp_volcano.scatter.json", x_title="Log2FC", y_title="D")
        else:
            self.chart_diffexp_volcano(cmps, ".diffexp", data, title, "diffexp_volcano.scatter.json")


        ctrl, test = cmps.split("_vs_")
        scatter_pd = diff_pd.loc[:, ['seq_id', ctrl + "_tpm", test + '_tpm', 'significant', 'regulate']]
        scatter_pd.rename(columns={ctrl + "_tpm": "group1", test + "_tpm": "group2"},inplace = True)


        # scatter_pd.set_index('seq_id', inplace=True)
        # scatter_pd = scatter_pd.loc[volcano_pd['seq_id'], :].reset_index()
        scatter_pd['group1'] = (scatter_pd['group1'] + 1).apply(np.log10)
        scatter_pd['group2'] = (scatter_pd['group2'] + 1).apply(np.log10)


        down_list_s = list()
        up_list_s = list()
        nosig_list_s = list()
        for vol in scatter_pd.to_dict('records'):
            if vol["regulate"] == "up" and vol["significant"] == "yes":
                up_list_s.append({
                    "y": vol["group2"],
                    "x": vol["group1"],
                    "selected": False,
                    "name": vol["seq_id"],
                    "color": "#FF2020"
                })
            elif vol["regulate"] == "down" and vol["significant"] == "yes":
                down_list_s.append({
                    "y": vol["group2"],
                    "x": vol["group1"],
                    "selected": False,
                    "name": vol["seq_id"],
                    "color": "#388E3C"
                })
            else:
                nosig_list_s.append({
                    "y": vol["group2"],
                    "x": vol["group1"],
                    "selected": False,
                    "name": vol["seq_id"],
                    "color": "#808080"
                })

        data = {
            "down": down_list_s,
            "up": up_list_s,
            "nosig": nosig_list_s
        }

        title = "{}.scatter".format(cmps)
        x_label = "log10({}_normalized + 1)".format(cmps.split("_vs_")[0])
        y_label = "log10({}_normalized + 1)".format(cmps.split("_vs_")[1])
        self.chart_diffexp_scatter2(cmps, ".diffexp", data, title, x_label, y_label, "diffexp_scatter.scatter.json")

    def chart_diffexp_volcano(self, name, out, data, title, json_mode, x_title=None, y_title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".volcano.js", 'w') as fo:
            a = json.loads(f.read())
            a["params"]["title"] = title
            a["data"] = data
            if x_title:
                a["params"]["x_label"] = x_title
            if y_title:
                a["params"]["y_label"] = y_title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".volcano.js", {"model": "highchart"}])

    def chart_diffexp_scatter2(self, name, out, data, title, x_label, y_label, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".scatter.js", 'w') as fo:
            a = json.loads(f.read())
            a["params"]["x_label"] = x_label
            a["params"]["y_label"] = y_label
            a["params"]["title"] = title
            a["data"] = data
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".scatter.js", {"model": "highchart"}])


    def chart_line_highchart(self, name, out, data, categories, json_mode, title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a["categories"] = categories
            a["data"] = data
            if title:
                a["params"]["title"] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {"model": "highchart", "highchart_type": "showCurve", "width": 650, "height": "430"}])

    def chart_columns_highchart(self, name, out, data, categories, json_mode, title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".columns.js", 'w') as fo:
            a = json.loads(f.read())
            a["categories"] = categories
            a["data"] = data
            if title:
                a["params"]["title"] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".columns.js", {"model": "highchart", "highchart_type": "showBar", "width": 650, "height": "430"}])

    def chart_snp_dis(self, snp_position_distribution, samples=None):
        snp_pos_pd =  pd.read_table(snp_position_distribution, header=0, sep='\t')
        headers = snp_pos_pd.columns
        samples = list()
        for sample in headers:
            if sample not in ['range_key', 'stat_type']:
                samples.append(sample)
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


    def chart_snp_stat(self, snp_stat, samples=None):
        snp_stat_pd =  pd.read_table(snp_stat, header=0, sep='\t')
        headers = snp_stat_pd.columns
        samples = list()
        for sample in headers:
            if sample not in ['range_key', 'stat_type']:
                samples.append(sample)
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
        headers = snp_stat_pd.columns
        samples = list()
        for sample in headers:
            if sample not in ['range_key', 'stat_type']:
                samples.append(sample)
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


    def chart_splice_all_stat(self, splice_all_stat, name="JC"):
        splice_pd =  pd.read_table(splice_all_stat, header=0, sep='\t')
        source = list()
        categories = list()
        for category in ["SE", "RI", "A5SS", "A3SS", "MXE"]:
            source.append({
                "category": category,
                "value": list(splice_pd[category])
            })

        categories = list(splice_pd["SAMPLE"])
        self.chart_column_stat( "all_" + name, "splice_stat", source, categories, "splice.all_stat.column.json")
        # 单样本饼图
        for row in splice_pd.iterrows():
            sample = row[1]["SAMPLE"]
            source =  [["value", "name", "category"]]
            for splice in ["MXE", "A5SS", "RI", "A3SS", "SE"]:
                if row[1]["TOTAL"] == 0:
                    continue
                else:
                    show = "{}:{}%".format(splice, "%0.4f" %(float(row[1][splice])/row[1]["TOTAL"] * 100))
                    source.append([row[1][splice], splice, show])

            title = "Summary of AS type"
            self.chart_pie("splice", "." + sample, source, "splice.stat.pie.json", title)

    def chart_column_stat(self, name, out, source, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".column.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".column.js", {}])

    def chart_column_color(self, name, out, source, color, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".column.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["series"][0]["visualMap"][0]["visualColorValue"] = color
            a["legend"][0]["color"] = color
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".column.js", {}])


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
            for category in ["SE", "RI", "A5SS", "A3SS", "MXE"]:
                sum_a += diff_stats_data[category + "_" + stat_type[stat]]
            for category in ["SE", "RI", "A5SS", "A3SS", "MXE"]:
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
                    "data": [psi_data[category + "_SAMPLE_1_" + f + "_exclusion"] for category in ["SE", "RI", "A5SS", "A3SS", "MXE"]],
                    "name": "exc"
                },
                {
                    "data": [psi_data[category + "_SAMPLE_1_" + f + "_inclusion"] for category in ["SE", "RI", "A5SS", "A3SS", "MXE"]],
                    "name": "inc"
                }
            ]
            title = "Summary of DES({})".format(cmp_name)
            self.chart_splice_diff_stat_column(cmp_name + "_" + stat, ".diff_splice_stat", source, splice_type, title, "splice.in_stat.column.json")

    def chart_splice_diff_stat_pie(self, name, out, source, categories, title=None, json_mode=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".pie.js", 'w') as fo:
            a = json.loads(f.read())
            if title:
                a["title"]["text"] = title
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".pie.js", {}])

    def chart_column(self, name, out, source, categories, title, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".column.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = title
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".column.js", {}])

    def chart_splice_diff_stat_column(self, name, out, source, categories, title, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".column.js", 'w') as fo:
            a = json.loads(f.read())
            if title:
                a["title"]["text"] = title
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".column.js", {}])

    def chart_class_column(self, name, out, source, categories, categories_source, title, json_mode, reset_margin=True):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".column.js", 'w') as fo:
            a = json.loads(f.read())
            if title:
                a["title"]["text"] = title
            a["dataset"][0]["source"] = source
            if categories:
                a["dataset"][0]["categories"] = categories
            if categories_source:
                a["dataset"][1]["source"] = categories_source

            if reset_margin:
                # 柱状图下边界加长
                if categories:
                    a = self.reset_margin(a, word_list = categories)
                elif categories_source:
                    a = self.reset_margin(a, word_list = source[0])
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".column.js", {"delay": 8000}])


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

        self.chart_column_stat("asprofile", ".stat", source, samples, "splice.all_stat.column.json")
        for rec in splice_pd.to_dict("records"):
            source = [["value", "name", "category"]]
            for category in categories:
                source.append([rec.get(category, 0),
                               category,
                               category + ":{}".format(rec.get(category))
                ])
            sample = rec["sample"]
            self.chart_pie("asprofile", "." + sample, source, "asprofile.stat.pie.json")

    def chart_gsva_diff(self, name, out, source, color, json_mode):
        json_mode = self.mode_dir + "/" + json_mode

        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".js", 'w') as fo:
            a = json.loads(f.read())
            a['series'][0]['visualMap'][0]['visualColorValue'] = color
            a['dataset'][0]['source'] = source
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".js", {}])

    def chart_stem(self, name, out, cluster_list, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]['source'] = cluster_list
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".js", {}])

    def chart_stem_profile(self, name, out, source_line, category_line, source_scatter, title, subtitle, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]['source'] = source_line
            a['dataset'][0]['categories'] = category_line
            a['dataset'][1]['source'] = source_scatter
            a["title"]["text"] = title
            a["title"]["title"]["text"] = title
            a['subtitle']['text'] = subtitle
            a["subtitle"]["subtitle"]["text"] = subtitle
            color_number = len(source_line)
            if color_number > 5:
                symbol = ["circle", "triangle", "diamond", "square", "triangle-down"]
                color = ["rgba(56,142,60,1)",
                        "rgba(244,67,54,1)",
                        "rgba(2,136,209,1)",
                        "rgba(255,152,0,1)",
                        "rgba(0,255,255,1)",
                        "rgba(233,30,99,1)",
                        "rgba(103,58,183,1)",
                        "rgba(0,100,0,1)",
                        "rgba(255,165,0,1)",
                        "rgba(255,0,0,1)",
                        "rgba(139,0,0,1)",
                        "rgba(128,128,128,1)",
                        "rgba(192,192,192,1)",
                        "rgba(144,238,144,1)",
                        "rgba(127,255,170,1)",
                        "rgba(255,20,147,1)",
                        "rgba(255,192,203,1)",
                        "rgba(255,0,255,1)",
                        "rgba(173,216,230,1)",
                        "rgba(135,206,235,1)", ]
                a["series"][0]["visualMap"][0]["visualSymbolValue"] = [symbol[i%5] for i in range(0, color_number)]
                a["series"][0]['visualMap'][0]['visualColorValue'] = [color[i%20] for i in range(0, color_number)]
                a["series"][1]["visualMap"][0]["visualSymbolValue"] = [symbol[i%5] for i in range(0, color_number)]
                a["series"][1]['visualMap'][0]['visualColorValue'] = [color[i%20] for i in range(0, color_number)]
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".js", {}])

    def _parse_dic_file(self, file):
        dic = dict(
            [(arr[0], int(arr[1])) for arr in [line.strip().split('\t') for line in open(file).readlines()[1:]]])
        return dic

    def chart_cog_class(self,geneset_list, cog_class_file_list):
        for geneset, cog_class_file in zip(geneset_list, cog_class_file_list):
            # print(geneset_list, cog_class_file_list)
            self.chart_geneset_class_cog(geneset, cog_class_file, geneset_list)
            # try:
            #     self.chart_geneset_class_cog(geneset,cog_class_file)
            # except:
            #     pass

    def chart_go_class(self,geneset_list, go_class_file_list):
        for geneset, go_class_file in zip(geneset_list, go_class_file_list):
            try:
                self.chart_diff_geneset_class_go(geneset,go_class_file)
            except:
                pass

    def chart_go_enrich(self,geneset_list, go_enrich_file_list):
        for geneset, go_enrich_file in zip(geneset_list, go_enrich_file_list):
            try:
                self.chart_diff_geneset_enrich_go(geneset,go_enrich_file,geneset_name = geneset)
            except:
                pass

    def chart_kegg_class(self,geneset_list,kegg_class_file_list,kegg_level_path):
        for geneset, kegg_class_file in zip(geneset_list, kegg_class_file_list):
            try:
                norm_kegg_class_file = self.convert_kegg_class_file(geneset,kegg_class_file,kegg_level_path)
                print norm_kegg_class_file
                self.chart_geneset_class_kegg(norm_kegg_class_file, [geneset])
            except Exception as e:
                print e
                pass

    def chart_kegg_enrich(self,geneset_list,kegg_enrich_file_list):
        for geneset, kegg_enrich_file in zip(geneset_list, kegg_enrich_file_list):
            try:
                self.chart_geneset_enrich_kegg(geneset,kegg_enrich_file,geneset_name = geneset)
            except:
                pass

    def chart_diff_genesets_venn(self, diff_genesets_outdir,genesets):
        # 通过获取差异基因集结果文件的原始ids文件来获取信息
        source_venn = list()
        prepare_json = json.load(open(os.path.join(diff_genesets_outdir, "prepare_json")))
        if len(prepare_json["genesets"]) <2:
            return

        for geneset, pre_dict in prepare_json["genesets"].items():
            geneset_dir = os.path.join(diff_genesets_outdir,geneset)
            geneset_list_file = pre_dict["geneset_path"] + "_gene.list"
            #modify by fwy 20210715 0个基因的基因集不在绘制venn图
            if int(pre_dict["gene_num"]) >0:
                with open(geneset_list_file,"r") as r:
                    genes = r.read().strip().split("\n")
                source_venn.append({
                    "data": genes,
                    "name": geneset
                })
        if len(source_venn) <=1:
            return
        self.chart_venn("diff_genesets", ".analysis", source_venn, "geneset.venn.venn.json")
        # self.chart_upset("diff_genesets", ".analysis", source_venn, "geneset.venn.upset.json")

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

        # 折线按表达量表排序
        if not samples_order:
            samples_order = list(cluster_pd.columns)

        cluster_pd["seq_id"] = cluster_pd.index

        seq_id = cluster_pd['seq_id'].tolist()
        for i in seq_id:
            if ':' in i:
                if i.replace(':', '-') in genes:
                    genes[genes.index(i.replace(':', '-'))] = i
            if '(' in i:
                if i.replace('(', '-').replace(')', '-') in genes:
                    genes[genes.index(i.replace('(', '-').replace(')', '-'))] = i

        cluster_pd = cluster_pd.loc[genes, :]

        corr_heat = [[""] + samples]
        for line_dict in cluster_pd.iterrows():
            corr_heat.append(
                [line_dict[0]] + [line_dict[1][s] for s in samples]
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
                             gene2group_source, "geneset.cluster.heatmap.json")

    def chart_geneset_class_cog(self, geneset_name, cog_class_table, geneset_list=None):
        a = pd.read_table(cog_class_table, header=0)
        col_names = a.columns
        # b = a.iloc[:, :-1]
        # b.columns = col_names[1:]
        b = a
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
        self.chart_class_column("{}.cog_annot".format(geneset_name), ".gene_set", source, categories, None, None, "geneset.annot_cog_stat.bar.json")

    def chart_diff_geneset_class_go(self,geneset_name, go_class_table, geneset_list=None, top=20):
        a = pd.read_table(go_class_table, header=0, index_col=0)
        # 原表格多一列处理
        col_names = a.columns
        b = a.iloc[:, :-1]
        b.columns = col_names[1:]
        b["sum"] = sum([b[x] for x in col_names if x.endswith("num")])

        b["k"] = range(b.shape[0])
        c = b.sort_values(by=["sum", "k"], ascending=[False, True])[:top]
        d = c.sort_index(ascending=False)
        d["type"] = d.index
        d = d.sort_values(by=["type", "sum", "k"], ascending=[False, False, True])
        categories = list(d["Term"])

        # c = b.sort_values(by=["sum"], ascending=[False])[:top]
        # d = c.sort_index(ascending=False)
        # d["type"] = d.index
        # d = d.sort_values(by=["type", "sum"], ascending=[False, False])
        # categories = list(d["Term"])
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
        #因为该交互工作流有且仅有一个基因集,所以直接取名字
        self.chart_class_column("{}.go_annot".format(geneset_name), ".gene_set", source, categories, class_source, None, bac_json, reset_margin=False)

    def chart_geneset_class_kegg(self, kegg_class_table, geneset_list=None, top=20):
        a = pd.read_table(kegg_class_table, header=0)
        # a = a[:top]

        geneset_list = [c.split("_genes")[0] for c in a.columns if c.endswith("genes")]
        for geneset in geneset_list:
            a_choose = a[a[geneset + "_genes"]>0]
            source = [
                ["item"] + list(a_choose["second_category"]),
                ["series"] + list(a_choose[geneset + "_genes"]),
                ["category"] + list(a_choose["first_category"])
            ]

            class_source = ["p1", "p2", "category"]
            for atype, sub_class in a.groupby("first_category"):
                class_source.append([sub_class.iloc[0,]["second_category"], sub_class.iloc[-1,]["second_category"], atype])
            print geneset
            self.chart_class_column("{}.kegg_annot.".format(geneset_list[0]),  geneset, source, None, class_source,  None,  "geneset.annot_kegg_stat.bar.json")



    def chart_diff_geneset_enrich_go(self,geneset, go_enrich_table, geneset_list=None, geneset_name="geneset",top=20,p_thre=0.5):
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
        self.chart_bar_and_line("{}.go_enrich".format(geneset), ".gene_set",  source_bar,  source_line, categories_line, geneset_name, "geneset.enrich_go.bar_line.json")

        b = a.sort_values(by=["go_type", "p_uncorrected"], ascending=[False, False])
        source_bar = [
            ["item"] + list(b["discription"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["go_type"])
        ]
        self.chart_bar("{}.go_enrich".format(geneset), ".gene_set",  source_bar, geneset_name, "geneset.enrich_go.bar.json")

        source_buble =  [["x", "y", "size", "fdr"]]
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["discription"],
                rec["study_count"],
                rec["p_corrected"],
            ])

        self.chart_buble("{}.go_enrich".format(geneset), ".gene_set",  source_buble, geneset_name, "geneset.enrich_go.buble.json")

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
        self.chart_buble2("{}.go_enrich".format(geneset), ".gene_set",  source_buble_list, geneset_name, "geneset.enrich_go.buble2.json")

    def chart_geneset_enrich_kegg(self, geneset,kegg_enrich_table, geneset_list=None, geneset_name="geneset", top=20,p_thre=0.5):
        a = pd.read_table(kegg_enrich_table, header=0)
        a["k"] = range(a.shape[0])
        a = a.sort_values(by=["Corrected P-Value", "k"], ascending=[True, True])[:top]
        a = a[a["Corrected P-Value"] <= p_thre]

        # a = a.sort_values(by=["P-Value"], ascending=[True])[:top]
        # a = a[a["Corrected P-Value"] <= p_thre]
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
        self.chart_bar_and_line("{}.kegg_enrich".format(geneset), ".gene_set",  source_bar,  source_line, categories_line, geneset_name, "geneset.enrich_kegg.bar_line.json")
        # kegg_order = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing', 'Cellular Processes', 'Organismal Systems', 'Human Diseases', 'Drug Development']

        #不用管这段代码干嘛的,有问题去看有参原版,这版为了适应医学的产品端zz需求弄的
        t= a[top-1::-1]
        kegg_order = list(t.drop_duplicates("typeI")["typeI"])
        def Reverse(lst):
            return [ele for ele in reversed(lst)]
        kegg_order = Reverse(kegg_order)
        t['order'] = [kegg_order.index(x) for x in t['typeI']]

        # kegg_order =list(a.drop_duplicates("typeI")["typeI"])
        # a['order']=[kegg_order.index(x) for x in a['typeI']]
        b = t.sort_values(by=["order", "neg_log10p_corrected"], ascending=[True, True])
        kegg_relate = {
            'Metabolism':"M",
            'Genetic Information Processing' :"GIP",
            'Environmental Information Processing' : "EIP",
            'Cellular Processes' : "CP",
            'Organismal Systems' : "OS",
            'Human Diseases' : "HD",
            "Drug Development" :"DD"
        }
        # b = a.sort_values(by=["P-Value"], ascending=[True])[:top]
        b["typeI"] = b.apply(lambda x: kegg_relate[x["typeI"]], axis=1)

        source_bar = [
            ["item"] + list(b["Term"]),
            ["series"] + list(b["neg_log10p_corrected"]),
            ["category"] + list(b["typeI"])
        ]
        self.chart_bar("{}.kegg_enrich".format(geneset), ".gene_set",  source_bar, geneset_name, "geneset.enrich_kegg.bar.json")

        source_buble =  [["x", "y", "size", "fdr"]]
        for rec in a.to_dict("records"):
            source_buble.append([
                rec["enrich_factor"],
                rec["Term"],
                rec["study_count"],
                rec["Corrected P-Value"]
            ])

        self.chart_buble("{}.kegg_enrich".format(geneset), ".gene_set",  source_buble, geneset_name, "geneset.enrich_kegg.buble.json")

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
        self.chart_buble2("{}.kegg_enrich".format(geneset), ".gene_set",  source_buble_list, geneset_name, "geneset.enrich_kegg.buble2.json")

    def convert_kegg_class_file(self,geneset,kegg_class_file,kegg_level_path):
        '''
        合并kegg通路统计与层级文件, 按level2统计基因数量
        '''
        if not os.path.exists(os.path.join("temporary",geneset)):
            os.makedirs(os.path.join("temporary",geneset))
        work_dir = os.path.join("temporary",geneset)
        path_stat_df = pd.read_table(kegg_class_file, header=0)
        path_level_df = pd.read_table(kegg_level_path, header=0)
        path_all = pd.merge(path_stat_df, path_level_df, on='Pathway_id')
        geneset_names = [x.split("_genes")[0] for x in path_all.columns if x.endswith("_genes")]


        # 按照kegg官网进行一级分类的排序
        list_custom = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing',
                       'Cellular Processes', 'Organismal Systems', 'Human Diseases',
                       'Drug Development']
        path_all['sort_k'] =  path_all['first_category'].map(lambda x:list_custom.index(x))
        path_all = path_all.sort_values(by="sort_k")
        path_all.drop(['graph_id', 'hyperlink', 'graph_png_id', 'sort_k'], axis=1, inplace=True)
        path_all.to_csv(work_dir + "/" + "kegg_annotation_analysis", sep='\t', index=False)

        # 按层级统计
        geneset_genes = [geneset_name + "_genes" for geneset_name in geneset_names]
        stat_df = path_all[['first_category', 'second_category'] + geneset_genes]
        stat_df.fillna("", inplace=True)
        def sum_gene(gene_str_list):
            gene_set = set()
            for genes_str in gene_str_list:
                for gene in genes_str.split(";"):
                    if gene== "":
                        continue
                    gene_set.add(gene.split("(")[0])
            return len(gene_set)

        agg_dict = {geneset_gene: sum_gene for geneset_gene in geneset_genes}
        stat_result = stat_df.groupby(['first_category', 'second_category']).agg(agg_dict)
        stat_result = stat_result.reset_index()
        stat_result['sort_k'] =  stat_result['first_category'].map(lambda x:list_custom.index(x))
        stat_result = stat_result.sort_values(by=["sort_k", 'second_category'])
        stat_result.drop(['sort_k'], axis=1, inplace=True)
        stat_result.to_csv(work_dir + "/" + "kegg_statistic", sep='\t', index=False)
        return work_dir + "/" + "kegg_statistic"

    #### 生成cmd脚本
    def generate_html_sh(self, para=True):
        self.command_list = list()
        for js in self.js_list:
            js_file = os.path.abspath(js[0])
            para_dict = js[1]
            html_file = os.path.splitext(js_file)[0] + ".html"
            pdf_file = os.path.splitext(js_file)[0] + ".pdf"
            if "model" in para_dict and para_dict["model"] == "highchart":
                mode_dir = self.mode_mode2
            else:
                mode_dir = self.mode_mode
            if "model" in para_dict and  para_dict["model"] == "highchart":
                if "use_medical" in para_dict and para_dict["use_medical"] == "yes":
                    html_mode = os.path.join(mode_dir, "sg_chart_model_medical.html")
                else:
                    html_mode = os.path.join(mode_dir, "sg_chart_model.html")
            else:
                if 'busco' in js_file:
                    html_mode = os.path.join(mode_dir, "sg_chart_model_busco.html")
                elif "use_new" in para_dict and para_dict["use_new"] == "yes":
                    html_mode = os.path.join(mode_dir, "sg_chart_model_new.html")
                else:
                    html_mode = os.path.join(mode_dir, "sg_chart_model.html")

            html = lxml.html.parse(html_mode)
            root = html.getroot()
            scripts = root.getchildren()
            for script in scripts:
                if script.tag == "body":
                    div1 = script.getchildren()[0]
                    # 修改图片打小
                    eles = div1.attrib['style'].split("; ")
                    if "height" in para_dict:
                        eles[0] = 'height: {}px'.format(para_dict["height"])
                    if "width" in para_dict:
                        eles[0] = 'width: {}px'.format(para_dict["width"])

                    div1.attrib['style'] = "; ".join(eles)

                if script.tag == "script" and "src" in script.attrib:
                    # print script.attrib
                    if script.attrib["src"] == "./sg_chart_model.js":
                        script.set('src', js_file)
                        script.text = ""
                    if script.attrib["src"] == "./sg.min.js":
                        if 'busco' in js_file:
                            script.set('src', os.path.join(mode_dir, "sg.min_tsg.js"))
                            script.text = ""
                        else:
                            pass
                    if script.attrib["src"].startswith("./"):
                        script.set('src', mode_dir + script.attrib["src"][1:])
                        script.text = ""


                else:
                    if "model" in para_dict and para_dict["model"] == "highchart":
                        if "highchart_type" in para_dict:
                            script.text = script.text.replace("showScatterMarkBig", para_dict["highchart_type"])
                        if "type" in para_dict and para_dict["type"] == "network":
                            script.text += "\n" + '$("#customize_button").remove();'
            html.write(html_file)

            if "delay" in para_dict:
                delay = para_dict["delay"]
            else:
                delay = 5000

            if "convert" in para_dict:
                cmd = "{}/wkhtmltopdf --no-stop-slow-scripts --enable-local-file-access --page-width {} --page-height {} --zoom {} {} {}".format(
                    self.wkhtml_dir,
                    para_dict["width"],
                    para_dict["height"],
                    para_dict["zoom"],
                    html_file,
                    pdf_file
                )
            elif "to_type" in para_dict and para_dict["to_type"] == "svg":
                svg_file = os.path.splitext(js_file)[0] + ".svg"
                cmd = "{}/phantomjs {}/sg_chart_phantome.svg.js {} {} && {} {} -o {}".format(
                    self.phantomjs_dir,
                    self.mode_mode,
                    html_file,
                    svg_file,
                    self.cairo_svg,
                    svg_file,
                    pdf_file
                )
            else:
                cmd = "{}/phantomjs {}/sg_chart_phantome.js {} {} {}".format(self.phantomjs_dir, self.mode_mode, html_file, pdf_file, delay)

            self.command_list.append(cmd)

        if para:
            with open(self.work_dir + "para_run.sh", 'w') as fo:
                fo.write("\n".join(self.command_list) + "\n")

    def para_to_pdf(self):
        '''
        # run in in tools
        '''
        pass

    def to_pdf(self):
        self.generate_html_sh(para=False)
        for command in self.command_list:
            print "command{}".format(command)
            os.system(command)
            # return_code = subprocess.call(command, shell=True)

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

        try:
            if "annot_stat" in chart_json:
                annot_stat = chart_json["annot_stat"]
                gene_exp = chart_json["gene_exp_all"]
                trans_exp = chart_json["tran_exp_all"]
                venn_dir = os.path.dirname(chart_json["annot_stat"])
                self.chart_annotation_stat("", gene_exp, trans_exp, annot_stat, venn_dir)
        except:
            pass
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
        '''
        # 差异图形在tools中生成
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

        if "snp_distribution" in chart_json:
            self.chart_snp_dis(chart_json["snp_distribution"])
        if "snp_stat" in chart_json:
            self.chart_snp_stat(chart_json["snp_stat"])
        if "snp_depth" in chart_json:
            self.chart_snp_depth_stat(chart_json["snp_depth"])

        #以下为差异一键化工作流
        if  "diff_geneset_venn" in chart_json:
            try:
                diff_genesets_out =  chart_json["diff_geneset_venn"]
                self.chart_diff_genesets_venn(diff_genesets_out,chart_json["genesets"])
            except Exception as e:
                print e
                pass

        if "cluster_geneset_name" in chart_json:
            self.chart_geneset_cluster(chart_json["cluster_exp"], chart_json["cluster_tree"], chart_json["sample_tree"],
                                       chart_json["subcluster_list"], group_dict=chart_json["group_dict"],seqid2namefile = chart_json["gene_annot_file"] )


        if "go_class" in chart_json:
            go_class_files = [chart_json["go_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_go_class(chart_json["genesets"], go_class_files)

        if "cog_class" in chart_json:
            cog_class_files = [chart_json["cog_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_cog_class(chart_json["genesets"], cog_class_files)

        if "go_enrich" in chart_json:
            go_enrich_files = [chart_json["go_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_go_enrich(chart_json["genesets"], go_enrich_files)

        if "kegg_class" in chart_json:
            # genesets = [chart_json["genesets"][0] + "_up", chart_json["genesets"][0] + "down"]
            genesets = chart_json["genesets"]
            kegg_class_files = [chart_json["kegg_class"].format(geneset_name=geneset) for geneset in genesets]

            self.chart_kegg_class(genesets, kegg_class_files, chart_json["kegg_level"])

        if "kegg_enrich" in chart_json:
            kegg_enrich_files = [chart_json["kegg_enrich"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
            self.chart_kegg_enrich(chart_json["genesets"], kegg_enrich_files)

        self.generate_html_sh()

if __name__ == '__main__':
    a = Chart()


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
    # a.chart_diffexp_scatter(diff_exp, "A_vs_B", soft="degseq")
    with open(sys.argv[1], 'r') as f:
        chart_json = json.load(f)
    genesets = chart_json["genesets"]
    kegg_class_files = [chart_json["kegg_class"].format(geneset_name=geneset) for geneset in genesets]

    a.chart_kegg_class(genesets, kegg_class_files, chart_json["kegg_level"])

    # genesets = chart_json["genesets"]
    # kegg_class_files = [chart_json["kegg_class"].format(geneset_name=geneset) for geneset in genesets]

    # self.chart_kegg_class(genesets, kegg_class_files, chart_json["kegg_level"])

    # diff_genesets_out =  chart_json["diff_geneset_venn"]
    # a.chart_diff_genesets_venn(diff_genesets_out,chart_json["genesets"])


    # if "cog_class" in chart_json:
    #     cog_class_files = [chart_json["cog_class"].format(geneset_name=geneset) for geneset in chart_json["genesets"]]
    #     a.chart_cog_class(chart_json["genesets"], cog_class_files)

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
