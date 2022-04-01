# -*- coding: utf-8 -*-
# __author__ = "zhangyitong"

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
from lxml import etree, html
from collections import OrderedDict
from collections import Counter
from bson.objectid import ObjectId
import glob
from pandas.api.types import CategoricalDtype


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
        self.mode_dir = chart_dir + "/prok_rna"
        self.mode_mode = chart_dir + "/model"
        self.mode_mode2 = chart_dir + "/prok_rna/highchart_model"
        self.js_list = list()
        self.work_dir = os.environ.get("PWD", "") + "/"
        try:
            os.link(self.mode_mode + '/iconfont.woff', self.work_dir + '/iconfont.woff')
        except:
            pass

    def prok_chart_raw_qc(self, sample_list, raw_qc_file_list, qc_type="raw"):
        for sample, qc_file in zip(sample_list, raw_qc_file_list):
            atcgn_pct = list()
            qc_file_list = qc_file.split(",")
            if len(qc_file_list) == 2:
                r_qc = pd.read_table(qc_file_list[0], header=0)
                r_qc['total_reads'] = r_qc["A_Count"] + r_qc["T_Count"] + r_qc["C_Count"] + r_qc["G_Count"] + r_qc["N_Count"]
                l_qc = pd.read_table(qc_file_list[1], header=0)
                l_qc['total_reads'] = l_qc["A_Count"] + l_qc["T_Count"] + l_qc["C_Count"] + l_qc["G_Count"] + l_qc["N_Count"]
                for base in ['A', 'T', 'C', 'G', 'N']:
                    r_qc["{}_pct".format(base)] = r_qc["{}_Count".format(base)]/r_qc["total_reads"] *100
                    l_qc["{}_pct".format(base)] = l_qc["{}_Count".format(base)] / l_qc["total_reads"] * 100
                r_qc["err"] = 10**(r_qc["mean"]/-10) * 100
                l_qc["err"] = 10**(l_qc["mean"]/-10) * 100

                source_base_data = [[i] + list(r_qc[i + "_pct"]) + list(l_qc[i + "_pct"]) for i in ["A", "T", "G", "C", "N"]]
                source_base_category = range(1, len(source_base_data[0]) - 1 + 1)

                source_err_data = [0.2 if i > 0.2 else i for i in list(r_qc["err"])] + [0.2 if i > 0.2 else i for i in list(l_qc["err"])]
                # source_err_data = list(r_qc["err"]) + list(l_qc["err"])
                source_err_data.insert(0, 'all')
                source_err_category = range(1, len(source_err_data))

                source_box_data_l = [[l_qc.iloc[i]["lW"],l_qc.iloc[i]["Q1"],l_qc.iloc[i]["med"],
                                      l_qc.iloc[i]["Q3"],l_qc.iloc[i]["rW"] ] for i in l_qc.index.tolist()]
                source_box_data_r = [[r_qc.iloc[i]["lW"],r_qc.iloc[i]["Q1"],r_qc.iloc[i]["med"],
                                      r_qc.iloc[i]["Q3"],r_qc.iloc[i]["rW"] ] for i in r_qc.index.tolist()]
                source_box_data = source_box_data_r + source_box_data_l
                source_box_category = range(1, len(source_box_data)+1)
            else:
                r_qc = pd.read_table(qc_file_list[0], header=0)
                r_qc['total_reads'] = r_qc["A_Count"] + r_qc["T_Count"] + r_qc["C_Count"] + r_qc["G_Count"] + r_qc[
                    "N_Count"]
                for base in ['A', 'T', 'C', 'G', 'N']:
                    r_qc["{}_pct".format(base)] = r_qc["{}_Count".format(base)]/r_qc["total_reads"] * 100
                r_qc["err"] = 10**(r_qc["mean"]/-10) * 100
                source_base_data = [[i] + list(r_qc[i + "_pct"]) for i in ["A", "T", "C", "G", "N"]]
                source_base_category = range(1, len(source_base_data[0]) - 1 + 1)
                source_err_data = [0.2 if i > 0.2 else i for i in list(r_qc["err"])]
                # source_err_data = list(r_qc["err"])
                source_err_data.insert(0, 'all')
                source_err_category = range(1, len(source_err_data))

                source_box_data_r = [[r_qc.iloc[i]["lW"],r_qc.iloc[i]["Q1"],r_qc.iloc[i]["med"],
                                      r_qc.iloc[i]["Q3"],r_qc.iloc[i]["rW"] ] for i in r_qc.index.tolist()]
                source_box_data = source_box_data_r
                source_box_category = range(1, len(source_box_data)+1)

            title_base = 'Nucleotide distribution of {} data:{}'.format(qc_type, sample)
            title_err = "Base error distribution of {} data:{}".format(qc_type, sample)
            title_qual = "Base quality distribution of {} data:{}".format(qc_type, sample)
            self.prok_chart_raw_qc_base(sample, ".{}_qc_base".format(qc_type), source_base_data, source_base_category, "qc_distribution_curve.json", title=title_base)
            self.prok_chart_raw_qc_err(sample, ".{}_qc_error".format(qc_type), source_err_data, source_err_category, "qc_error_rate.json", title=title_err)
            self.prok_chart_raw_qc_qual(sample, ".{}_qc_qual".format(qc_type), source_box_data, source_box_category, "qc_mass_distribution.json", title=title_qual)

    def prok_chart_raw_qc_base(self, name, out, data_list, category_list, json_mode, title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'] = data_list
            a['categories'] = category_list
            if title:
                a["params"]["title"] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {"model": "highchart", "highchart_type": "showCurve", "width": 500, "height": "400"}])

    def prok_chart_raw_qc_err(self, name, out, data_list, category_list, json_mode, title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'] = [data_list]
            a['categories'] = category_list
            if title:
                a["params"]["title"] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {"model": "highchart", "highchart_type": "showCurve", "width": 500, "height": "400"}])

    def prok_chart_raw_qc_qual(self, name, out, data_list, category_list, json_mode, title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".box.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'] = data_list
            a['categories'] = category_list
            if title:
                a["params"]["title"] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".box.js", {"model": "highchart", "highchart_type": "showBoxPlot", "width": 500, "height": "400"}])

    def prok_map_saturation(self, sample_list, saturate_files):
        for sample, s_file in zip(sample_list, saturate_files):
            r_file = s_file + '.saturation.R'
            percent_file = s_file + '.cluster_percent.xls'
            data_list = list()
            with open(r_file, "r") as f:
                for line in f:
                    if re.match(r"legend", line):
                        all_num = re.findall("num=[\d]*", line)
                        categories = [
                            "[0-0.3)=" + all_num[0][4:],
                            "[0.3-0.6)=" + all_num[1][4:],
                            "[0.6-3.5)=" + all_num[2][4:],
                            "[3.5-15)=" + all_num[3][4:],
                            "[15-60)=" + all_num[4][4:],
                            ">=60=" + all_num[5][4:],
                        ]
            with open(percent_file, 'r') as p:
                for line in p:
                    items = line.strip().split()[1:]
                    data_list.append([float(x) for x in items])
            title = 'Saturation curve of sequencing({})'.format(sample)
            self.prok_chart_map_assessment(sample, ".map_saturation", data_list, categories, "assessment_saturation_curve.json", title=title)

    def prok_map_coverage(self, coverage_files):
        sample_list = list()
        data_list = list()
        for each in coverage_files:
            sample_name = os.path.basename(each).split(".")[0][9:]
            sample_list.append(sample_name)
            with open(each, "r") as f:
                f.readline()
                value = f.next().strip().split()
                plot_value = list()
                for i in range(100):
                    plot_value.append(int(eval(value[i + 1])))
            data_list.append(plot_value)
        colours = ['#388E3C' for i in sample_list]
        self.prok_chart_map_assessment('', 'map_coverage', data_list, sample_list, 'assessment_coverage.json', colour=colours)

    def prok_chart_map_assessment(self, name, out, data_list, category_list, json_mode, title=None, colour=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'] = data_list
            a['legend'] = category_list
            if title:
                a["params"]["title"] = title
            if colour:
                a['params']['colors'] = colour
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {"model": "highchart", "highchart_type": "showCurve", "width": 650, "height": "410"}])

    def prok_annot_venn(self, venn_dir):
        data_list_venn = []
        venn_db = {
            'nr': 'NR',
            'swissprot': 'Swiss-Prot',
            'pfam': 'Pfam',
            'kegg': 'KEGG',
            'go': 'GO',
            'cog': 'COG',
            'subloc': 'SubCell-Location',
        }
        db_order = [i+'_venn.txt' for i in ['nr', 'swissprot', 'pfam', 'cog', 'go', 'kegg', 'subloc']]
        for each in sorted(os.listdir(venn_dir), key=lambda x: db_order.index(x)):
            db_name = each.split('_venn')[0]
            if db_name in venn_db.keys():
                with open(os.path.join(venn_dir, each), 'r') as f:
                    gene_list = f.read().strip().split('\n')
                data_list_venn.append({"data": gene_list, "name": venn_db[db_name]})
        self.prok_chart_highchart_venn('', 'stats_annot_venn', data_list_venn, 'ref_annot_venn.json')

    def prok_annot_bar(self, stats_file):
        stats_pd = pd.read_table(stats_file, header=0)
        # categories = stats_pd[stats_pd.columns[0]].values.tolist()[:-2]
        # data_list_num = stats_pd[stats_pd.columns[1]].values.tolist()[:-2]
        categories = ['NR', 'Swiss-Prot', 'Pfam', 'COG', 'GO', 'KEGG'] # ordered
        data_list_num = [stats_pd[stats_pd[stats_pd.columns[0]] == i].values[0][1] for i in categories]
        data_list_percent = [(stats_pd[stats_pd[stats_pd.columns[0]] == i].values[0][2])*100 for i in categories]
        self.prok_chart_highchart_bar('', 'stats_annot_num_bar', data_list_num, 'ref_annot_bar.json', categories=categories)
        self.prok_chart_highchart_bar('', 'stats_annot_percent_bar', data_list_percent, 'ref_annot_bar.json',
                                      categories=categories, y_lab='Percent of Unigenes (%)')

    def prok_chart_highchart_venn(self, name, out, data_list, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".venn.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'] = data_list
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".venn.js",
                                 {"model": "highchart", "highchart_type": "showVenn", "width": 650, "height": "430", "type": "annot_venn"}])

    def prok_chart_highchart_bar(self, name, out, data_list, json_mode, categories=None, x_lab=None, y_lab=None, legend=None, dimensions=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".bar.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'][0] = data_list
            if categories:
                a['categories'] = categories
            if x_lab:
                a['params']['x_label'] = x_lab
            if y_lab:
                a['params']['y_label'] = y_lab
            if legend:
                a['legend'] = legend
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            if dimensions:
                width = dimensions['width']
                height = dimensions['height']
            else:
                width = 650
                height = '430'
            self.js_list.append([self.work_dir + name + out + ".bar.js",
                                 {"model": "highchart", "highchart_type": "showBar", "width": width, "height": height}])

    def prok_annot_cog_bar(self, cog_file):
        data_list = list()
        categories = list()
        cog_df = pd.read_table(cog_file, header=0)
        cog_list = cog_df.iloc[:, [2, 3]].values.tolist()
        for each in cog_list:
            data_list.append(each[1])
            categories.append(each[0])
        self.prok_chart_annot_cog_bar('', 'annot_cog', data_list, categories, 'ref_cog_class_bar.json')

    def prok_chart_annot_cog_bar(self, name, out, data_list, categories, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".cog_bar.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'][0] = data_list
            a['categories'] = categories
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".cog_bar.js",
                                 {"model": "highchart", "highchart_type": "showBar", "width": 650, "height": "430"}])

    def prok_annot_go_bar_pie(self, level2_file):
        level2_df = pd.read_table(level2_file, header=0, dtype={'Percent': 'str'})
        bar_bp_top10 = level2_df[level2_df[level2_df.columns[0]] == 'biological_process'].iloc[:10, [0, 2, 3, 4]].values.tolist()   #GO (Lev1), GO Term, Seq Number, Percent
        bar_cc_top10 = level2_df[level2_df[level2_df.columns[0]] == 'cellular_component'].iloc[:10,[0, 2, 3, 4]].values.tolist()
        bar_mf_top10 = level2_df[level2_df[level2_df.columns[0]] == 'molecular_function'].iloc[:10,[0, 2, 3, 4]].values.tolist()
        data_list_bar = bar_bp_top10 + bar_cc_top10 + bar_mf_top10

        pie_bp_top10 = [{'name': '{} : {}({})'.format(index+1, item[1], str(item[2])), 'value': item[2]} for index, item in enumerate(bar_bp_top10)]  # {"name": "9 : zinc ion binding(74)","value": 74}
        pie_cc_top10 = [{'name': '{} : {}({})'.format(index+1, item[1], str(item[2])), 'value': item[2]} for index, item in enumerate(bar_cc_top10)]
        pie_mf_top10 = [{'name': '{} : {}({})'.format(index+1, item[1], str(item[2])), 'value': item[2]} for index, item in enumerate(bar_mf_top10)]
        data_list_pie = [pie_bp_top10, pie_cc_top10, pie_mf_top10]

        self.prok_chart_annot_go_bar('', 'annot_go_bar', data_list_bar, 'ref_go_class_bar.json')
        self.prok_chart_annot_go_pie('', 'annot_go_pie', data_list_pie, 'ref_go_class_pie.json')

    def prok_chart_annot_go_bar(self, name, out, data_list, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".multi_bar.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'] = data_list
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4, sort_keys=True))
            self.js_list.append([self.work_dir + name + out + ".multi_bar.js",
                                 {"model": "highchart", "highchart_type": "go_bar", "width": 1000, "height": "430"}])

    def prok_chart_annot_go_pie(self, name, out, data_list, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".multi_pie.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'] = data_list
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".multi_pie.js",
                                 {"model": "highchart", "highchart_type": "multiPie", "width": 1000, "height": "430"}])

    def prok_annot_kegg_bar(self, kegg_file=None, kegg_table=None, name=''):
        if kegg_file:
            kegg_df = pd.read_table(kegg_file, header=None)
        else:
            kegg_df = kegg_table
        cat_order = CategoricalDtype(['Metabolism', 'Genetic Information Processing',
                                      'Environmental Information Processing', 'Cellular Processes',
                                      'Organismal Systems', 'Human Diseases', 'Drug Development'], ordered=True)
        kegg_df[kegg_df.columns[0]] = kegg_df[kegg_df.columns[0]].astype(cat_order)
        kegg_df.sort_values(by=[kegg_df.columns[0], kegg_df.columns[1]], inplace=True)
        category_list = ['category'] + kegg_df[kegg_df.columns[0]].tolist()
        des_list = ['item'] + kegg_df[kegg_df.columns[1]].tolist()
        num_list = ['series'] + kegg_df[kegg_df.columns[2]].tolist()
        source_list = [des_list, num_list, category_list]
        summary_dict = OrderedDict()
        for index, item in enumerate(category_list[1:]):
            if item not in summary_dict.keys():
                summary_dict[item] = list()
            summary_dict[item].append(index)
        summary_list = [["p1", "p2", "category"]] + [[min(value), max(value), key] for key, value in summary_dict.items() if value]
        self.prok_chart_annot_kegg_bar(name, 'annot_kegg_bar', source_list, summary_list, 'ref_kegg_pathway.json')

    def prok_chart_annot_kegg_bar(self, name, out, data_list, summary_list, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".kegg_bar.js", 'w') as fo:
            a = json.loads(f.read())
            a['dataset'][0]['source'] = data_list
            a['dataset'][1]['source'] = summary_list
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4, sort_keys=True))
            self.js_list.append([self.work_dir + name + out + ".kegg_bar.js", {"width": 900, "height": "500"}])

    def process_exp_matrix(self, exp_matrix, log_base=None, group_dict=None):
        if type(exp_matrix) == str or type(exp_matrix) == bytes or isinstance(exp_matrix, unicode):
            all_exp_pd = pd.read_table(exp_matrix, index_col=0, header=0)
        else:
            print(exp_matrix, 'is assumed to be a pandas DataFrame Object')
            all_exp_pd = exp_matrix
        all_exp_pd.index.name = 'seq_id'
        if 'type' in all_exp_pd.columns.tolist():
            temp_type = all_exp_pd.loc[:, 'type']
            tmp_genename = all_exp_pd.loc[:, 'gene_name']
            all_exp_pd.drop(['type', 'gene_name', 'biotype'], axis=1, inplace=True)

        if group_dict is not None:
            group_exp = list()
            for g in group_dict:
                g_exp = all_exp_pd.loc[:, group_dict[g]].mean(axis=1)
                g_exp.name = g
                group_exp.append(g_exp)
            all_exp_pd = pd.concat(group_exp, axis=1)

        if log_base:
            if 'type' in all_exp_pd.columns.tolist():
                all_exp_pd = all_exp_pd.applymap(lambda x: float(x))
                if log_base == math.e:
                    all_exp_pd = np.log(all_exp_pd + 1)
                    all_exp_pd['type'] = temp_type
                    all_exp_pd['gene_name'] = tmp_genename
                elif log_base == 2:
                    all_exp_pd = np.log2(all_exp_pd + 1)
                    all_exp_pd['type'] = temp_type
                    all_exp_pd['gene_name'] = tmp_genename
                elif log_base == 10:
                    all_exp_pd = np.log10(all_exp_pd + 1)
                    all_exp_pd['type'] = temp_type
                    all_exp_pd['gene_name'] = tmp_genename
                else:
                    raise Exception('log base of {} is not supported'.format(log_base))
            else:
                all_exp_pd = all_exp_pd.applymap(lambda x: float(x))
                if log_base == math.e:
                    all_exp_pd = np.log(all_exp_pd + 1)
                elif log_base == 2:
                    all_exp_pd = np.log2(all_exp_pd + 1)
                elif log_base == 10:
                    all_exp_pd = np.log10(all_exp_pd + 1)
                else:
                    raise Exception('log base of {} is not supported'.format(log_base))

        if not group_dict and not log_base:
            all_exp_pd['type'] = temp_type
            all_exp_pd['gene_name'] = tmp_genename
        return all_exp_pd

    def get_box_source(self, exp_pd, order):
        # 获取box source
        stat_list = self.get_box(exp_pd)
        source_all = []
        for each in order:
            for box in stat_list:
                if box['sample'] == each:
                    source_all.append({
                        "category": box['sample'],
                        "data": [
                            box['min'],
                            box['q1'],
                            box['median'],
                            box['q3'],
                            box['max']
                        ],
                        "name": box['sample']
                    })
        return source_all

    def get_box(self, all_exp_pd):
        """
        get box plot info for each column of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        """
        stat_dict_list = list()
        if 'type' in all_exp_pd.columns.tolist():
            temp_type = all_exp_pd.loc[:, 'type']
            tmp_genename = all_exp_pd.loc[:, 'gene_name']
            all_exp_pd.drop(['type', 'gene_name'], axis=1, inplace=True)
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
                'min-q1': lt25,
                'q1-median': lt50 - lt25,
                'median-q3': lt75 - lt50,
                'q3-max': exp_pd.shape[0] - lt75,
                'upper_whisker': upper_whisker,
                'lower_whisker': lower_whisker,
                'upper_outliers': upper_outliers,
                'lower_outliers': lower_outliers,
            })
            stat_dict_list.append(tmp_dict)
        return stat_dict_list

    def get_density(self, all_exp_pd):
        """
        sampling 1000 density point for each columns of the input pandas DataFrame
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        """
        records = list()
        if 'type' in all_exp_pd.columns.tolist():
            temp_type = all_exp_pd.loc[:, 'type']
            tmp_genename = all_exp_pd.loc[:, 'gene_name']
            all_exp_pd.drop(['type', 'gene_name'], axis=1, inplace=True)

        target_columns = all_exp_pd.columns
        for sample in target_columns:
            exp = all_exp_pd[sample]
            exp = exp[exp != 0]
            if len(exp) > 1:
                density_func = stats.gaussian_kde(exp)
                min_exp, max_exp = exp.min(), exp.max()
                x_data = np.linspace(min_exp, max_exp, num=1000, endpoint=False)
                y_data = density_func(x_data)
                point_dict_list = pd.DataFrame({'log10exp': x_data, 'density': y_data}) # 修改这个字段的名字，因为实际是取的10的对数
                records.append(dict(sample=sample, data=point_dict_list))
            else:
                min_exp, max_exp = 0, 0
                x_data = np.linspace(min_exp, max_exp, num=1000, endpoint=False)
                y_data = density_func(x_data)
                point_dict_list = pd.DataFrame({'log10exp': x_data, 'density': y_data})
                records.append(dict(sample=sample, data=point_dict_list))
        return records

    def get_density_source(self, exp_pd, order):
        # 获取box source
        stat_list = self.get_density(exp_pd)
        source_all = []
        for each in order:
            for density in stat_list:
                if density['sample'] == each:
                    source_all.append({
                        "data": list(density["data"]["density"]), # 有点问题每个样本的横坐标不是一个体系，category如何设置
                        "name": density['sample']
                    })

        category = list(stat_list[0]["data"]["log10exp"])
        return category, source_all

    def get_violin_source(self, exp_pd, order):
        # 获取box source
        source_all = [["data", "name"]]
        target_columns = exp_pd.columns[1:]
        for sample in order:
            exp = exp_pd[sample]
            source_all.append([list(exp), sample])
        return source_all

    def chart_exp_dis(self, gene_exp, group_dict, samples):
        # 使用参考转录本做
        # if os.path.exists(tran_exp):
        #     levels = ["G", "T"]
        # else:
        for classify in ["samples", "groups"]:
            if classify == "samples":
                g = None
                s = samples
                order = samples
            else:
                g = group_dict
                s = None
                order = group_dict.keys()
            exp = gene_exp

            all_exp_pd = self.process_exp_matrix(exp, log_base=10, group_dict=g)

            # all_exp_pd = all_exp_pd[all_exp_pd.sum(axis=1) > 0.001]
            source_all = self.get_box_source(all_exp_pd, order)
            self.chart_raw_exp_box(classify, ".exp_distribution", source_all, "exp_distribution_box.json")

            category, source_density = self.get_density_source(all_exp_pd, order)
            self.chart_raw_exp_density(classify, ".exp_distribution", source_density, category, "exp_distribution_density.json")

            all_exp_pd = all_exp_pd.sample(frac=0.7)
            all_exp_pd.reset_index(level=0, inplace=True)
            source_violin = self.get_violin_source(all_exp_pd, order)
            self.chart_raw_exp_violin(classify, ".exp_distribution", source_violin, "exp_distribution_violin.json", y_title=None)

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

        self.chart_exp_venn_venn("all", ".exp", source_venn, "exp_all_venn.json", bar_show=True)

    def chart_exp_venn_venn(self, name, out, source, json_mode, bar_show=False):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".venn.js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]["source"] = source

            if bar_show:
                a['series'][0]['bar_show'] = bar_show
            geneset_list = [s['name'][10:] for s in source]
            a = self.reset_margin(j=a, margin_type="right", word_list=geneset_list)
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".venn.js",
                                 {"convert": "wkhtmltopdf", "zoom": 7, "width": 800, "height": 600}])

    def chart_exp_corr(self, sample_corr_file, sample_tree, group_dict=None):
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
            group_dict = json.loads(group_dict)
            for g, ss in group_dict.items():
                for s in ss:
                    sample2group_source.append([s, g])
        else:
            for s in samples:
                sample2group_source.append([s, s])

        self.chart_heat_tree("exp", ".heatmap", corr_source, sample_tree, sample_tree, sample2group_source, sample2group_source, "exp_all_heatmap.json")

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

        self.chart_exp_pca_scatter("all", ".exp_relation_pca", pca_source, x_title, y_title, "exp_all_pca.json")

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

            self.chart_exp_pca_scatter2("all", ".exp_relation_pca_ell", pca_source, ellipse_source, x_title, y_title, "exp_all_pca2.json")

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
            if len(categories) > 2:
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
            a['series'][1]['visualMap'][0]['visualColorValue'] = ['white' for i in source_ellipse]
            categories = set([s[-1] for s in source])
            if len(categories) > 2:
                symbol = ["circle", "triangle", "diamond", "square", "triangle-down"]
                a["series"][0]["visualMap"][0]["visualSymbolValue"] = [symbol[i%5] for i in range(0, len(categories))]
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".scatter.js", {}])

    def chart_diffexp_stat(self, summary, cmp_list=None, soft=None):
        with open(summary, "r") as f:
            lines = f.readlines()
            if len(lines) < 3:
                return

        summary_pd = pd.read_table(summary, header=[0, 1])
        levels = summary_pd.columns.levels
        labels = summary_pd.columns.labels
        summary_pd.columns = levels[0][labels[0]]

        if cmp_list:
            categories = ['_vs_'.join(i) for i in cmp_list]
        else:
            categories = sorted(list(summary_pd.columns[1:-1]))
        up_list = list()
        down_list = list()
        for cmp1 in categories:
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
        self.chart_diffexp_stat_bar(soft, ".diffexp_summary", all_source, categories, "diff_exp_stats_bar.json")
        self.chart_diffexp_stat_bar2(soft, ".diffexp_summary", all_source2, categories, "diff_exp_stats_stackedbar.json")

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
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".stacked_bar.js", 'w') as fo:
            a = json.loads(f.read())
            a["title"]["text"] = a["title"]["text"].format(sample_name = name)
            a["dataset"][0]["source"] = source
            a["dataset"][0]["categories"] = categories
            a = self.reset_margin(j=a, margin_type="bottom", word_list=categories)
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".stacked_bar.js", {}])

    def prok_diffexp_scatter(self, diff_exp, cmps, soft=None, pvalue_padjust='padjust'):
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
            self.prok_chart_diffexp_volcano(cmps, ".diffexp", data, title, "diff_exp_volcano.json", x_title="Log2FC", y_title="D")
        else:
            self.prok_chart_diffexp_volcano(cmps, ".diffexp", data, title, "diff_exp_volcano.json")

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
        x_label = "log10({}_mean) TPM".format(cmps.split("_vs_")[0])
        y_label = "log10({}_mean) TPM".format(cmps.split("_vs_")[1])
        self.prok_chart_diffexp_scatter2(cmps, ".diffexp", data, title, x_label, y_label, "diff_exp_scatter.json")

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

    def prok_chart_diffexp_volcano(self, name, out, data, title, json_mode, x_title=None, y_title=None):
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
            self.js_list.append([self.work_dir + name + out + ".volcano.js", {"model": "highchart", "highchart_type": "showScatterMarkBig", "width": 650, "height": "430"}])

    def prok_chart_diffexp_scatter2(self, name, out, data, title, x_label, y_label, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".scatter.js", 'w') as fo:
            a = json.loads(f.read())
            a["params"]["x_label"] = x_label
            a["params"]["y_label"] = y_label
            a["params"]["title"] = title
            a["data"] = data
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".scatter.js", {"model": "highchart", "highchart_type": "showScatterMarkBig", "width": 650, "height": "430"}])

    def chart_geneset_venn_ids(self, geneset_ids):
        project_type = 'prok_rna'
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

        self.chart_exp_venn_venn("geneset", "",  source, "geneset_venn.json", bar_show=True)

    def chart_geneset_venn_file(self, geneset_file):
        all_detail = pd.read_table(geneset_file, header=0, sep='\t')
        all_detail['type'] = np.where(all_detail['seq_id'].str.contains("sRNA"), 'sRNA', 'mRNA')
        all_detail.fillna('', inplace=True)
        geneset_num = len(set(all_detail['compare']))
        if geneset_num > 6:
            cmp_except = list(set(all_detail['compare']))[6:]
        else:
            cmp_except = list()
        venn_source_list = list()

        for i in all_detail.groupby('compare'):
            compare = i[0]
            if cmp_except and compare in cmp_except:
                continue
            diff_pd = i[1]
            name = compare.replace('|', '_vs_') + '_mRNA'
            sig_seqs = diff_pd[(~diff_pd['seq_id'].str.contains("sRNA")) & (diff_pd['significant'] == 'yes')]['seq_id'].tolist()
            venn_source_list.append({"data": sig_seqs, "name": name})
            if geneset_num < 3:
                sig_seqs_up = diff_pd[(~diff_pd['seq_id'].str.contains("sRNA")) & (diff_pd['significant'] == 'yes') & (diff_pd['regulate'] == 'up')]['seq_id'].tolist()
                sig_seqs_down = diff_pd[(~diff_pd['seq_id'].str.contains("sRNA")) & (diff_pd['significant'] == 'yes') & (diff_pd['regulate'] == 'down')]['seq_id'].tolist()
                venn_source_list.append({"data": sig_seqs_up, "name": name + '_up'})
                venn_source_list.append({"data": sig_seqs_down, "name": name + '_down'})

        self.chart_exp_venn_venn("geneset", "",  venn_source_list, "geneset_venn.json", bar_show=True)

    def chart_geneset_venn(self, geneset_ids):
        print(geneset_ids)
        project_type = 'prok_rna'
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
        self.chart_exp_venn_venn("geneset", "", source, "geneset_venn.json", bar_show=True)

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

            sub_pd = cluster_pd.loc[sub_genes, :]
            mean_df = pd.DataFrame([[sub_pd[col].mean() if col != 'seq_id' else 'mean' for col in sub_pd.columns]],
                                   columns=sub_pd.columns, index=['mean'])      # work out mean in each sample
            sub_pd = sub_pd.append(mean_df)
            sub_data_list = list()
            legend_list = list()
            for rec in sub_pd.to_dict("records"):
                # print rec
                sub_data_list.append([rec[s] for s in samples_order])
                legend_list.append(rec["seq_id"])
            # print source_line[:3]
            title = "subcluster_{}({} genes)".format(c_num, len(sub_pd)-1)
            colours = ['#F0F0F0' for i in range(0, len(sub_pd)-1)] + ['#000099']
            categories = samples_order
            self.chart_line_point("subcluster", "." + c_num, sub_data_list, legend_list, categories, title, colours, "geneset_subcluster_line.json")

        for gene in genes:
            if gene in gene2group_dict:
                gene2group_source.append([gene, gene2group_dict[gene]])

        self.chart_heat_tree("geneset", ".cluster", corr_heat, sample_tree, gene_tree,  sample2group_source, gene2group_source, "geneset_cluster_heatmap.json")

    def chart_line_point(self, name, out, source, source_point, categories, title, colour, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".line.js", 'w') as fo:
            a = json.loads(f.read())
            a["params"]["title"] = title
            a["data"] = source
            a["categories"] = categories
            a["legend"] = source_point
            a['params']['colors'] = colour
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".line.js", {"model": "highchart", "highchart_type": "showCurve", 'height': 400, 'width': 700}])

    def chart_geneset_class_cog(self, cog_class_table, geneset_list=None):
        a = pd.read_table(cog_class_table, header=0, index_col=0)
        col_names = a.columns
        b = a.iloc[:, :-1]
        b.columns = col_names[1:]
        b.sort_values(by='Functional Categoris', ascending=True, inplace=True)
        categories = [x.split()[0][1] for x in list(b["Functional Categoris"])]
        if categories[0] == ':':
            categories = [x.split(':')[0] for x in list(b["Functional Categoris"])]
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
        self.chart_class_column("cog_annot", ".gene_set", source, categories, None, None, "geneset_cog_bar.json")

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

    def prok_geneset_goclass(self, geneset_goclass_file):
        geneset_go_df = pd.read_table(geneset_goclass_file, index_col=False)
        # ____________单向直方图-up_down__________
        colnum = int(geneset_go_df.shape[1]) / 3
        need_list_1 = [(i + 1) * 3 - 3 for i in range(colnum)][1:]  # [3,6]
        need_list_2 = list(
            set(list(range(colnum * 3))) - set([(i + 1) * 3 - 1 for i in range(colnum)]))  # [0,1,3,4,6,7]
        need_list_3 = [(i + 1) * 2 + 1 for i in range(colnum - 1)]  # [3,5]
        print need_list_1, need_list_2, need_list_3

        sub_cat = [list(geneset_go_df.columns)[i].split(' ')[0] for i in need_list_1]
        for i in range(len(sub_cat)):
            sort_cols = [0, 1, (i + 1) * 3, (i + 1) * 3 + 1]
            geneset_go_df_top20 = geneset_go_df.sort_values(
                by=[geneset_go_df.columns[need_list_1[i]], geneset_go_df.columns[1]],
                ascending=[False, True]).iloc[:, sort_cols].drop_duplicates(subset=None, keep='first',
                                                                            inplace=False).iloc[:20, :]
            geneset_go_df_top20[geneset_go_df_top20.columns[3]] = \
            geneset_go_df_top20[geneset_go_df_top20.columns[3]].map(
                lambda x: x.split('(')[0]).astype(float)
            geneset_go_df_top20 = geneset_go_df_top20.sort_values(by=[geneset_go_df_top20.columns[0], geneset_go_df_top20.columns[2], geneset_go_df_top20.columns[1]], ascending=[True, False, True])
            list_list = geneset_go_df_top20.values.tolist()
            self.prok_chart_geneset_goclass(sub_cat[i] + '.go_bar_geneset', list_list, "geneset_go_class_bar.json", y_right=True)

        if len(sub_cat) > 1:
            geneset_go_df['sum'] = 0
            for j in need_list_1:
                geneset_go_df['sum'] += geneset_go_df.iloc[:, j]
            geneset_go_df_top20 = geneset_go_df.sort_values(by=['sum', geneset_go_df.columns[1]], ascending=[False, True]).iloc[:, need_list_2].drop_duplicates(subset=None, keep='first',inplace=False).iloc[:20, :]
            for ii in need_list_3:
                geneset_go_df_top20[geneset_go_df_top20.columns[ii]] = \
                geneset_go_df_top20[geneset_go_df_top20.columns[ii]].map(
                    lambda x: x.split('(')[0]).astype(float)
            list_list = geneset_go_df_top20.values.tolist()
            self.prok_chart_geneset_goclass('go_bar_geneset', list_list, "geneset_go_class_bar.json", sub_cat=sub_cat)

    def prok_chart_geneset_goclass(self, name, list_list, json_mode, sub_cat=None, y_right=False):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + ".go_bar.js", 'w') as fo:
            a = json.loads(f.read())
            a["data"] = list_list
            if sub_cat:
                a["sub_categories"] = sub_cat
            if y_right:
                a['params']["yAxis2Display"] = True
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4, sort_keys=True))
            self.js_list.append([self.work_dir + name + ".go_bar.js",
                                 {"model": "highchart", "highchart_type": "go_bar", "width": 1100, "height": 600}])

    def prok_geneset_kegg_bar(self, kegg_file, kegg_level):
        level_df = pd.read_table(kegg_level, header=0)
        level_pair = level_df[['Pathway_id', 'second_category', 'first_category']]
        geneset_df = pd.read_table(kegg_file, header=0)
        geneset_info = [[0, 1, i, i+1] for i in range(2, len(geneset_df.columns)-1, 2)]     # [i,i+1]=[[2,3],[4,5],etc.]
        for each in geneset_info:
            each_df = geneset_df[[geneset_df.columns[i] for i in each]]
            each_df = each_df[each_df[each_df.columns[2]] != 0]
            if each_df.empty:
                continue
            geneset_name = each_df.columns[2].split('_numbers')[0]
            geneset_kegg = level_pair.merge(each_df, left_on='Pathway_id', right_on=each_df.columns[0])
            geneset_col = geneset_kegg.columns[-1]
            geneset_uniq = geneset_kegg.groupby('second_category')[geneset_col].apply(lambda x: x.str.cat(sep=';')).reset_index()
            level_info = geneset_kegg[['second_category', 'first_category']].drop_duplicates()
            geneset_filtered = pd.DataFrame(geneset_uniq).merge(level_info, on='second_category')
            geneset_filtered['count'] = geneset_filtered.apply(lambda x: len(set([each.split('(')[0] for each in x[geneset_col].split(';')])), axis=1)
            geneset_filtered.sort_values(['first_category'], inplace=True)
            geneset_table = geneset_filtered.iloc[:, [2, 0, 3]]
            self.prok_annot_kegg_bar(kegg_table=geneset_table, name=geneset_name + '.')

    def chart_geneset_gorich(self, go_enrich_file, geneset_name=''):
        # for colour bar
        colour_bar = OrderedDict()
        colour_list = ['#FF2020', '#80FB39', '#35F3F8', '#7503EB']
        for i, j in enumerate(['0%', '40%', '65%', '100%']):
            colour_bar[j] = colour_list[i]

        # parse file
        go_enrich_all_file_df = pd.read_table(go_enrich_file)
        if 'p_corrected' not in go_enrich_all_file_df.columns:
            go_enrich_all_file_df['p_corrected'] = self.multtest_correct(go_enrich_all_file_df['p_uncorrected'].tolist(), "bh")

        # Bar pdf
        go_enrich_all_file_df = go_enrich_all_file_df[(go_enrich_all_file_df.iloc[:, 2] == 'e')]
        if not go_enrich_all_file_df.empty:
            go_enrich_all_file_df['new_percent'] = go_enrich_all_file_df.iloc[:, 4].str.split(pat="/", expand=True).iloc[:, 0].astype(int) / go_enrich_all_file_df.iloc[:, 5].str.split(pat="/", expand=True).iloc[:, 0].astype(int)
            go_enrich_all_data_ = go_enrich_all_file_df.sort_values(by=[go_enrich_all_file_df.columns[-2], go_enrich_all_file_df.columns[3]], ascending=True).iloc[:20, [1, 3, -1, -2]]
            go_enrich_all_data = go_enrich_all_data_.values.tolist()
            self.chart_geneset_gorichBar(geneset_name + '_go_enrich_bar', go_enrich_all_data, "geneset_enrich_go_shadowbar.json", geneset_name, colour_bar)

        ## DenseBubble pdf
        def get_sort(x, b_type='scatter'):
            df = x.sort_values(by=x.columns[3], ascending=True)
            if b_type == 'dense':
                data_list = [[i[0], eval(i[1].split('/')[0]), eval(i[2].split('/')[0]), i[3]] for i in
                             df.values.tolist()]
            else:
                data_list = [{'y': i[5], 'x': i[-1], 'desc': i[0], 'name': i[6], 'size': eval(i[1].split('/')[0])} for i
                             in df.values.tolist()]
            return data_list

        def bubble_each(go_enrich_file_df, data_type, group_name):
            pvalues = go_enrich_file_df['p_uncorrected'].values.tolist()
            if len([i for i in pvalues if i > 0]) > 0:
                pvalues_min = min([i for i in pvalues if i > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = -math.log10(pvalues_min)
            go_enrich_file_df['logUnpvalue'] = go_enrich_file_df.apply(lambda x: -math.log10(x['p_uncorrected']) if x['p_uncorrected'] > 0 else pvalues_min, axis=1)
            go_enrich_dense_data_ = go_enrich_file_df.sort_values(by=[go_enrich_file_df.columns[-3], go_enrich_file_df.columns[1]], ascending=True).iloc[:20, [3, 4, 5, -3, 1, -1, 0, -2]]    # name, ratio_in_study, ratio_in_pop, p_corrected, NS, logUnpvalue, GO, new_percent
            go_enrich_dense_data = get_sort(go_enrich_dense_data_, 'dense')

            go_enrich_bubble_data_ = go_enrich_file_df.sort_values(by=['p_uncorrected', go_enrich_file_df.columns[1]], ascending=True).iloc[:20,[3, 4, 5, -3, 1, -1, 0, -2]]  # name, ratio_in_study, ratio_in_pop, p_corrected, NS, logUnpvalue, GO, new_percent
            go_enrich_all_data__ = go_enrich_bubble_data_
            bp_data = get_sort(go_enrich_all_data__[go_enrich_all_data__[go_enrich_all_data__.columns[4]] == 'BP'])
            cc_data = get_sort(go_enrich_all_data__[go_enrich_all_data__[go_enrich_all_data__.columns[4]] == 'CC'])
            mf_data = get_sort(go_enrich_all_data__[go_enrich_all_data__[go_enrich_all_data__.columns[4]] == 'MF'])
            go_enrich_scatter_data = [bp_data, cc_data, mf_data]

            self.chart_geneset_gorichDenseBubble(group_name + '_'+data_type + '_go_enrich_dense_bubble', go_enrich_dense_data, False, 'densebubble', "geneset_enrich_go_densebubble.json", group_name, colour_bar=colour_bar)
            self.chart_geneset_gorichDenseBubble(group_name + '_'+data_type + '_go_enrich_scatter_bubble', go_enrich_scatter_data, True, 'scatterbubble', "geneset_enrich_go_scatterbubble.json", group_name)
        if not go_enrich_all_file_df.empty:
            bubble_each(go_enrich_all_file_df, 'all', geneset_name)

    def chart_geneset_gorichBar(self, name, data, json_mode, group_name, colour_bar):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + ".shadowbar.js", 'w') as fo:
            a = json.loads(f.read())
            a["data"] = data
            a["params"]["text"] = a["params"]["text"].format(group_name)
            a['params']['linearGradient'] = colour_bar
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + ".shadowbar.js",
                                 {"model": "highchart", "highchart_type": "shadow_bar", "width": 850, "height": 650}])

    def chart_geneset_gorichDenseBubble(self, name, data, scatter_or_not, bubble_type, json_mode, group_name, colour_bar=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + "." + bubble_type + ".js", 'w') as fo:
            a = json.loads(f.read())
            a["data"] = data
            a['params']['disperse'] = scatter_or_not
            a["params"]["title"] = a["params"]["title"].format(group_name)
            if colour_bar:
                a['params']['linearGradient'] = colour_bar
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + "." + bubble_type + ".js",
                                 {"model": "highchart", "highchart_type": "graph.bubble", "width": 750, "height": 550}])

    def chart_geneset_keggrich(self, kegg_enrich_all_file, geneset_name):

        abbrev_dic = {'Cellular Processes': 'CP', 'Environmental Information Processing': 'EIP',
                      'Genetic Information Processing': 'GIP', 'Human Diseases': 'HD', 'Metabolism': 'M',
                      'Organismal Systems': 'OS'}

        # for colour bar
        colour_bar = OrderedDict()
        colour_list = ['#FF2020', '#80FB39', '#35F3F8', '#7503EB']
        for i, j in enumerate(['0%', '40%', '65%', '100%']):
            colour_bar[j] = colour_list[i]

        # Bar pdf
        kegg_enrich_all_file_df = pd.read_table(kegg_enrich_all_file)
        kegg_enrich_all_file_df['new_percent'] = kegg_enrich_all_file_df.iloc[:, 4].str.split(pat="/",
                                                                                              expand=True).iloc[:,
                                                 0].astype(int) / kegg_enrich_all_file_df.iloc[:, 5].str.split(
            pat="/", expand=True).iloc[:, 0].astype(int)
        kegg_enrich_all_data_ = kegg_enrich_all_file_df.sort_values(
            by=[kegg_enrich_all_file_df.columns[7], kegg_enrich_all_file_df.columns[11]],
            ascending=[True, False]).iloc[:20, [11, 1, 12, 7]]
        kegg_enrich_all_data_ = kegg_enrich_all_data_[kegg_enrich_all_data_['Corrected P-Value'] <= 0.5]
        kegg_enrich_all_data_['typeI'] = kegg_enrich_all_data_['typeI'].map(abbrev_dic)
        kegg_enrich_all_data = kegg_enrich_all_data_.values.tolist()
        self.chart_geneset_keggrichBar(geneset_name + '_enrichkegg', kegg_enrich_all_data, "geneset_enrich_kegg_shadowbar.json", geneset_name, colour_bar)

        # Dense Bubble
        kegg_enrich_all_data_ = kegg_enrich_all_file_df.sort_values(by=[kegg_enrich_all_file_df.columns[7], kegg_enrich_all_file_df.columns[11]], ascending=True).iloc[:20,[1, 4, 5, 7]]
        kegg_enrich_all_data_ = kegg_enrich_all_data_[kegg_enrich_all_data_['Corrected P-Value'] <= 0.5]
        kegg_enrich_all_data = [[i[0], eval(i[1].split('/')[0]), eval(i[2].split('/')[0]), i[3]] for i in kegg_enrich_all_data_.values.tolist()]

        self.chart_geneset_keggrichDenseBubble(geneset_name + '_enrichkegg', kegg_enrich_all_data, False, 'densebubble',
                                               "geneset_enrich_kegg_densebubble.json",geneset_name, colour_bar)

    def chart_geneset_keggrichBar(self, name, data, json_mode, geneset_name, colour_bar):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + ".shadowbar.js", 'w') as fo:
            a = json.loads(f.read())
            a["data"] = data
            a['params']['linearGradient'] = OrderedDict(
                (k, a['params']['linearGradient'].get(k)) for k in ["0%", "40%", "65%", "100%"])
            a['params']["text"] = a['params']["text"].format(geneset_name)
            a['params']['linearGradient'] = colour_bar
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + ".shadowbar.js",
                                 {"model": "highchart", "highchart_type": "shadow_bar", "width": 850, "height": 550}])

    def chart_geneset_keggrichDenseBubble(self, name, data, scatter_or_not, bubble_type, json_mode, geneset_name, colour_bar):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + "." + bubble_type + ".js", 'w') as fo:
            a = json.loads(f.read())
            a["data"] = data
            a['params']['disperse'] = scatter_or_not
            a['params']["title"] = a['params']["title"].format(geneset_name)
            a['params']['linearGradient'] = colour_bar
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + "." + bubble_type + ".js",
                                 {"model": "highchart", "highchart_type": "graph.bubble", "width": 900, "height": 550}])

    def prok_geneset_circ(self, geneset_circ_choose, geneset_circ_zscore, annot_type):
        if os.path.exists(geneset_circ_choose) and os.path.exists(geneset_circ_zscore):
            circ_pd = pd.read_table(geneset_circ_choose)
            body = circ_pd.iloc[:, :-1].values.tolist()
            pdf_height_index = pd.read_table(geneset_circ_choose).shape[0]
            header = pd.read_table(geneset_circ_zscore, header=None, index_col=0).values.tolist()
            self.prok_chart_geneset_circ('chord', body, header, pdf_height_index, "geneset_enrich_circ.json", annot_type)

    def prok_chart_geneset_circ(self, name, body, header, pdf_height_index, json_mode, annot_type):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + ".circ.js", 'w') as fo:
            a = json.loads(f.read())
            a["table"]["body"] = body
            a["table"]["header"] = header
            if annot_type == 'kegg':
                a['title'] = 'KEGG Pathway'
                a['legendTitle'] = 'KEGG Pathway'
            else:
                a['title'] = 'GO Term'
                a['legendTitle'] = 'GO Term'
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            # 41*2*k=838*3.14, k=32, so height=(rows*2*32)/3.14=rows*20.38~~rows*21
            self.js_list.append([self.work_dir + name + ".circ.js",
                                 {"model": "highchart", "highchart_type": "circ", "width": 65 * pdf_height_index,
                                  "height": 40 * pdf_height_index}])

    def prok_srna_predict_length(self, predict_file):
        predict_df = pd.read_table(predict_file, header=0)
        length_list = predict_df['length'].values.tolist()
        step = 50
        categories = range(0, 501, step)
        predict_list = [[j for j in length_list if categories[i] < j <= categories[i+1]] for i in range(len(categories))[:-1]]
        predict_list.append([j for j in length_list if categories[-1] < j])
        data_list = [len(i) for i in predict_list]
        self.prok_chart_highchart_bar('', 'srna_length_bar', data_list, 'srna_length_bar.json', dimensions={'width': 1000, 'height': '600'})

    def prok_srna_stats_venn(self, list_dir):
        list_files = [glob.glob(os.path.join(list_dir, '*_vs_{}.list'.format(i)))[0] for i in ['sRNATarBase', 'SIPHI', 'sRNAMap', 'rfam', 'BSRD']]
        data_list = list()
        for file in list_files:
            annot_type = os.path.basename(file).strip('.list').split('_vs_')[-1]
            if annot_type.lower() == 'rfam':
                annot_type = 'Rfam'
            with open(file, 'r') as f:
                srna_list = [i for i in f.read().split('\n') if i]
            if annot_type.lower() != 'bsrd' and len(srna_list) > 0:
                data_list.append({'data': srna_list, "name": annot_type})
        self.prok_chart_highchart_venn('', 'srna_annot_venn', data_list, 'srna_annot_venn.json')

    def prok_srna_rfam_pie(self, rfam_file):
        rfam_df = pd.read_table(rfam_file, header=None)
        rfam_df.drop(rfam_df[(rfam_df[0] == 'All_matched') | (rfam_df[0] == 'All_reads')].index, inplace=True)
        data_list = [{'y': i[1], 'name':i[0]} for i in rfam_df.values.tolist()]
        self.prok_chart_highchart_pie('', 'srna_rfam_pie', data_list, 'srna_annot_rfam_pie.json')

    def prok_chart_highchart_pie(self, name, out, data_list, json_mode, title=None):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".pie.js", 'w') as fo:
            a = json.loads(f.read())
            a['data'] = data_list
            if title:
                a['params']['title'] = title
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4))
            self.js_list.append([self.work_dir + name + out + ".pie.js",
                                 {"model": "highchart", "highchart_type": "showPie", "width": 650, "height": "430"}])

    def prok_operon_bar(self, operon_file):
        operon_pd = pd.read_table(operon_file, header=0)

        # frequency of operon length
        operon_pd['len'] = operon_pd['Stop'] - operon_pd['Start'] + 1
        len_list = operon_pd['len'].values.tolist()
        step = 2000
        categories = self.get_categories(max(len_list), step)
        length_list = [[j for j in len_list if categories[i] < j <= categories[i + 1]] for i in range(len(categories))[:-1]]
        data_list = [len(i) for i in length_list]
        cat_list = [str(categories[i]+1) + '~' + str(categories[i+1]) for i in range(len(categories))[:-1]]
        self.prok_chart_highchart_bar('', 'operon_len_bar', data_list, 'operon_length_bar.json', categories=cat_list, dimensions={'width': 1000, 'height': '700'})

        # distribution of gene number
        num_list = operon_pd['Number_of_genes'].values.tolist()
        num_categories = range(0, max(num_list) + 1)    # step=1
        gene_num_list = [num_list.count(i) for i in num_categories[1:]]
        num_cat = [str(i) for i in num_categories[1:]]
        self.prok_chart_highchart_bar('', 'operon_gene_num_bar', gene_num_list, 'operon_gene_num_bar.json', categories=num_cat, dimensions={'width': 1000, 'height': '700'})

    def get_categories(self, max_num, step, limit=None):
        if limit and max_num > limit:
            categories = range(0, limit + 1, step)
        elif int(max_num/step) > 0 and max_num % step > 0:
            categories = range(0, (int(max_num/step) + 1) * step + 1, step)
        elif int(max_num/step) == 0:
            categories = [0, step]
        else:
            categories = range(0, max_num + 1, step)
        return categories

    def prok_utr_len_bar(self, utr_file):

        def get_data(utr_len, utr_step):
            utr_small = [i for i in utr_len if i <= 1000]
            utr_large = [i for i in utr_len if i > 1000]
            utr_categories = self.get_categories(max(utr_small), utr_step)
            utr_list = [[j for j in utr_small if utr_categories[i] < j <= utr_categories[i + 1]] for i in range(len(utr_categories))[:-1]]
            cat_utr = [str(utr_categories[i] + 1) + '~' + str(utr_categories[i + 1]) for i in range(len(utr_categories))[:-1]]
            if utr_large:
                utr_list.append(utr_large)
                cat_utr.append('>1000')
            data_utr = [len(i) for i in utr_list]
            return data_utr, cat_utr

        utr_pd = pd.read_table(utr_file, header=0)
        utr_pd['len'] = abs(utr_pd['End'] - utr_pd['Start']) + 1
        step = 50

        # bar chart for Distribution of 5'UTR length
        utr5_len = utr_pd[utr_pd['Type'] == 'UTR5']['len'].values.tolist()
        if utr5_len:
            data_utr5, cat_utr5 = get_data(utr5_len, step)
            self.prok_chart_highchart_bar('UTR5_', 'len_bar', data_utr5, 'utr_length_bar.json', categories=cat_utr5, dimensions={'width': 700, 'height': 700})

        # bar chart for Distribution of 3'UTR length
        utr3_len = utr_pd[utr_pd['Type'] == 'UTR3']['len'].values.tolist()
        if utr3_len:
            data_utr3, cat_utr3 = get_data(utr3_len, step)
            self.prok_chart_highchart_bar('UTR3_', 'len_bar', data_utr3, 'utr_length_bar.json', categories=cat_utr3, dimensions={'width': 700, 'height': 700})

    def prok_snp_pie_bar(self, mutation_file, sample_list):
        mut_df = pd.read_table(mutation_file, header=0)
        mut_df['snp_type'] = mut_df['REF'] + "/" + mut_df['ALT']
        mut_df['type'] = mut_df.apply(lambda x: "snp" if len(x['snp_type']) == 3 and "-" not in x['snp_type'] else "indel", axis=1)

        anno_type = set(mut_df['ANNO'])
        title = '{} distribution in genome regions ({})'
        title_pie = 'Summary of SNPs ({})'
        snp_df = mut_df[mut_df['type'] == 'snp']
        for sample in sample_list:
            # SNP/InDel distribution of genomic regions
            sample_df = mut_df[~mut_df[sample].isin(['./.', '0/0'])]
            snp_anno = sample_df[sample_df['type'] == 'snp']['ANNO'].values.tolist()
            indel_anno = sample_df[sample_df['type'] == 'indel']['ANNO'].values.tolist()
            snp_list = [{'y': snp_anno.count(i), 'name': i} for i in anno_type]
            indel_list = [{'y': indel_anno.count(i), 'name': i} for i in anno_type]
            self.prok_chart_highchart_pie(sample, '.snp_regions_pie', snp_list, 'snp_indel_distribution_pie.json', title=title.format('SNP', sample))
            self.prok_chart_highchart_pie(sample, '.indel_regions_pie', indel_list, 'snp_indel_distribution_pie.json',
                                          title=title.format('InDel', sample))

            # depth stats
            depth_list = mut_df.apply(lambda x: int(x[sample].split("/")[1]), axis=1).tolist()
            depth_categories = self.get_categories(max(depth_list), 100, limit=500)
            depth_categories.insert(1, 30)
            depth_stats_list = [[j for j in depth_list if depth_categories[i] < j <= depth_categories[i + 1]] for i in
                                range(len(depth_categories))[:-1]]
            depth_stats_list.append([j for j in depth_list if depth_categories[-1] < j])
            depth_list_bar = [len(i) for i in depth_stats_list]
            cat_depth = ['<=30'] + [str(depth_categories[i] + 1) + '-' + str(depth_categories[i + 1]) for i in
                                    range(len(depth_categories))[1:-1]] + ['>500']
            depth_list_pie = [{'y': i, 'name': j} for i, j in zip(depth_list_bar, cat_depth)]
            self.prok_chart_highchart_bar(sample, '.depth_stats_bar', depth_list_bar, 'snp_stats_bar.json',
                                          categories=cat_depth, x_lab='Depth', y_lab='SNP Number', legend=[sample])
            self.prok_chart_highchart_pie(sample, '.depth_stats_pie', depth_list_pie, 'snp_stats_pie.json', title=title_pie.format(sample))

            # frequency stats
            frequency_list = sample_df[sample_df['ANNO'] == 'exonic'].iloc[:, 7].tolist()
            frequency_dict = Counter(frequency_list)
            frequency_count = dict(frequency_dict).values()
            cat_freq = [1, 2, 3, 4, '>=5']
            freq_list = [frequency_count.count(i) for i in cat_freq[:-1]]
            freq_list.append(len([j for j in frequency_count if 4 < j]))
            freq_list_pie = [{'y': i, 'name': j} for i, j in zip(freq_list, cat_freq)]
            self.prok_chart_highchart_bar(sample, '.frequency_stats_bar', freq_list, 'snp_stats_bar.json',
                                          categories=cat_freq, x_lab='SNP number per gene', y_lab='Gene Number', legend=[sample])
            self.prok_chart_highchart_pie(sample, '.frequency_stats_pie', freq_list_pie, 'snp_stats_pie.json',
                                          title=title_pie.format(sample))

            # type stats
            cat_type = ["A/G", "A/C", "C/T", "G/A", "G/C", "C/A", "A/T", "C/G", "G/T", "T/C", "T/A", "T/G"]
            type_list = sample_df[sample_df['type'] == 'snp']['snp_type'].tolist()
            type_list_bar = [type_list.count(i) for i in cat_type]
            type_list_pie = [{'y': i, 'name': j} for i, j in zip(type_list_bar, cat_type)]
            self.prok_chart_highchart_bar(sample, '.type_stats_bar', type_list_bar, 'snp_stats_bar.json',
                                          categories=cat_type, x_lab='SNP type', y_lab='Number of SNPs',
                                          legend=[sample])
            self.prok_chart_highchart_pie(sample, '.type_stats_pie', type_list_pie, 'snp_stats_pie.json',
                                          title=title_pie.format(sample))

    def chart_stem_cluster(self, profile_list, cluster_table, method):
        if method == 'SCM':
            cluster = pd.read_table(cluster_table, header=0, index_col=None, sep='\t', dtype=object)

            cluster_new = cluster[['profile', 'model', 'cluster', 'pvalue']]
            model_dict = OrderedDict(cluster_new[['profile', 'model']].values.tolist())
            cluster_new['pvalue'] = cluster_new['pvalue'].astype('float')   # sorting data by pvalue
            cluster_new.sort_values(by=['pvalue'], inplace=True)
            cluster_new['pvalue'] = cluster_new['pvalue'].astype(object)
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
                # title = "Profile #{}(0.0,-1.0,-1.0,0.0,0.0,-1.0)".format(profile_num)
                title = "Profile #{}({})".format(profile_num, model_dict[str(profile_num)])
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

    def chart_stem(self, name, out, cluster_list, json_mode):
        json_mode = self.mode_dir + "/" + json_mode
        with open(json_mode, 'r') as f, open(self.work_dir + name + out + ".js", 'w') as fo:
            a = json.loads(f.read())
            a["dataset"][0]['source'] = cluster_list
            fo.write("var options = ")
            fo.write(json.dumps(a, indent=4, sort_keys=True))
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

    def multtest_correct(self, p_values, methods='bonferroni'): #复制于package下的kegg_rich.py
        """
        1. Bonferroni
        2. Bonferroni Step-down(Holm)
        3. Benjamini and Hochberg False Discovery Rate
        4. FDR Benjamini-Yekutieli
        :param pvalue_list:
        :param methods:
        :return: np.array
        """

        def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):
            '''pvalue correction for false discovery rate
            This covers Benjamini/Hochberg for independent or positively correlated and
            Benjamini/Yekutieli for general or negatively correlated tests. Both are
            available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.
            Parameters
            ----------
            pvals : array_like
                set of p-values of the individual tests.
            alpha : float
                error rate
            method : {'indep', 'negcorr')
            Returns
            -------
            rejected : array, bool
                True if a hypothesis is rejected, False if not
            pvalue-corrected : array
                pvalues adjusted for multiple hypothesis testing to limit FDR
            Notes
            -----
            If there is prior information on the fraction of true hypothesis, then alpha
            should be set to alpha * m/m_0 where m is the number of tests,
            given by the p-values, and m_0 is an estimate of the true hypothesis.
            (see Benjamini, Krieger and Yekuteli)
            The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
            of false hypotheses will be available (soon).
            Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
            fdr_by.
            '''

            def _ecdf(x):
                '''
                no frills empirical cdf used in fdrcorrection
                '''
                nobs = len(x)
                return np.arange(1, nobs + 1) / float(nobs)

            pvals = np.asarray(pvals)
            if not is_sorted:
                pvals_sortind = np.argsort(pvals)
                pvals_sorted = np.take(pvals, pvals_sortind)
            else:
                pvals_sorted = pvals  # alias

            if method in ['i', 'indep', 'p', 'poscorr']:
                ecdffactor = _ecdf(pvals_sorted)
            elif method in ['n', 'negcorr']:
                cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))  # corrected this
                ecdffactor = _ecdf(pvals_sorted) / cm
            else:
                raise ValueError('only indep and negcorr implemented')
            reject = pvals_sorted <= ecdffactor * alpha
            if reject.any():
                rejectmax = max(np.nonzero(reject)[0])
                reject[:rejectmax] = True

            pvals_corrected_raw = pvals_sorted / ecdffactor
            pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
            del pvals_corrected_raw
            pvals_corrected[pvals_corrected > 1] = 1
            if not is_sorted:
                pvals_corrected_ = np.empty_like(pvals_corrected)
                pvals_corrected_[pvals_sortind] = pvals_corrected
                del pvals_corrected
                reject_ = np.empty_like(reject)
                reject_[pvals_sortind] = reject
                return reject_, pvals_corrected_
            else:
                return reject, pvals_corrected

        pvalue_list = list(p_values)
        n = len(pvalue_list)
        fdr = list()
        if methods == 'bonferroni':
            fdr = [eachP * n for eachP in pvalue_list]
        elif methods == 'holm':
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP * (n - sorted_pvalues.index(eachP)) for eachP in pvalue_list]
        elif methods == 'bh':
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP * n / (sorted_pvalues.index(eachP) + 1) for eachP in pvalue_list]
        elif methods == 'by':
            _, fdr = fdrcorrection(pvalue_list, alpha=0.05, method='negcorr', is_sorted=False)
        fdr = np.array(fdr)
        fdr[fdr > 1] = 1.
        return fdr

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
            if "Histogram_of_KEGG" in js[0]:
                html_mode = os.path.join(mode_dir, "sg_chart_model_drna_KEGG.html")
            # elif "violin" in js[0]:
            #     html_mode = os.path.join(mode_dir, "sg_chart_model_drna_violin.html")
            # elif "exp.venn" in js[0] or "exp.upset" in js[0] or "exp_relation_pca" in js[0] or "exp.heat_corr" in js[0]:
            # elif 'all.exp' in js[0] or 'violin' in js[0]:
            #     html_mode = os.path.join(mode_dir, "sg_chart_model.html")
            elif "annot_venn" in js[0]:
                html_mode = os.path.join(mode_dir, "sg_chart_model_jvenn.html")
            elif 'go_bar' in js[0]:
                html_mode = os.path.join(mode_dir, 'sg_chart_model_go_bar.html')
            elif 'go_pie' in js[0]:
                html_mode = os.path.join(mode_dir, 'sg_chart_model_multiPie.html')
            elif 'shadowbar' in js[0]:
                html_mode = os.path.join(mode_dir, 'sg_chart_model_shadow_bar.html')
            elif 'densebubble' in js[0]:
                html_mode = os.path.join(mode_dir, 'sg_chart_model_bubble_dense.html')
            elif 'scatterbubble' in js[0]:
                html_mode = os.path.join(mode_dir, 'sg_chart_model_bubble_scatter.html')
            elif 'chord.circ' in js[0]:
                html_mode = os.path.join(mode_dir, 'sg_chart_model_enrich_circ.html')
            else:
                html_mode = os.path.join(mode_dir, "sg_chart_model.html")
            html = lxml.html.parse(html_mode)
            root = html.getroot()
            scripts = root.getchildren()
            for script in scripts:
                if script.tag == "body":
                    div1 = script.getchildren()[0]
                    # 修改图片打小
                    if "model" in para_dict and para_dict["model"] == "highchart":
                        eles = div1.attrib['style'].strip(' ').split(';')
                        if "height" in para_dict:
                            eles.append('height: {}px'.format(para_dict["height"]))
                        if "width" in para_dict:
                            eles.append('width: {}px'.format(str(para_dict["width"])))
                    else:
                        eles = list()
                        if "width" in para_dict:
                            eles.append('min-width: {}px'.format(str(para_dict["width"])))
                        if "height" in para_dict:
                            eles.append('height: {}px'.format(para_dict["height"]))
                    div1.attrib['style'] = ";".join([i for i in eles if i])

                if script.tag == "script" and "src" in script.attrib:
                    # print script.attrib
                    if script.attrib["src"] == "./sg_chart_model.js":
                        script.set('src', js_file)
                        script.text = ""
                        # if "Histogram_of_KEGG" in js[0]:
                        #     script.set('src', os.path.join(mode_dir, 'bar_line.js'))
                        #     script.text = ""
                    elif script.attrib["src"].startswith("./"):
                        script.set('src', mode_dir + script.attrib["src"][1:])
                        script.text = ""
                else:
                    if "model" in para_dict and para_dict["model"] == "highchart":
                        if "highchart_type" in para_dict:
                            if para_dict['highchart_type'] != 'go_bar':
                                script.text = script.text.replace("showScatterMarkBig", para_dict["highchart_type"])
                            else:
                                pass
                        if "type" in para_dict and para_dict["type"] in ["network"]:
                            script.text += "\n" + '$("#customize_button").remove();'

            # if "Histogram_of_KEGG" in js[0]:
            #     script_new = etree.Element('script')
            #     script_new.set('src', os.path.join(mode_dir, 'bar_line.js'))
            #     script_new.text = ""
            #     scripts += [script_new]

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
            else:
                cmd = "{}/phantomjs {}/sg_chart_phantome.js {} {} {}".format(self.phantomjs_dir, self.mode_mode, html_file, pdf_file, delay)

            self.command_list.append(cmd)

        if para:
            with open(self.work_dir + "para_run.sh", 'w') as fo:
                fo.write("\n".join(self.command_list) + '\n')

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
        if "qc_file_before" in chart_json:
            qc_files = [chart_json["qc_file_before"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.prok_chart_raw_qc(chart_json["samples"], qc_files)
        if "qc_file_after" in chart_json:
            qc_files = [chart_json["qc_file_after"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.prok_chart_raw_qc(chart_json["samples"], qc_files, qc_type="clean")
        if 'map_saturation' in chart_json:
            sat_files = [chart_json["map_saturation"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.prok_map_saturation(chart_json["samples"], sat_files)
        if 'map_coverage' in chart_json:
            cov_files = [chart_json["map_coverage"].format(sample_name=sample) for sample in chart_json["samples"]]
            self.prok_map_coverage(cov_files)
        if 'annot_venn' in chart_json:
            self.prok_annot_venn(chart_json['annot_venn'])
        if 'annot_stats' in chart_json:
            self.prok_annot_bar(chart_json['annot_stats'])
        if 'cog_annot' in chart_json:
            self.prok_annot_cog_bar(chart_json['cog_annot'])
        if 'go_annot_level2' in chart_json:
            self.prok_annot_go_bar_pie(chart_json['go_annot_level2'])
        if 'kegg_annot' in chart_json:
            self.prok_annot_kegg_bar(kegg_file=chart_json['kegg_annot'])
        if 'exp_matrix' and 'group_dict' in chart_json:
            self.chart_exp_dis(chart_json['exp_matrix'], chart_json['group_dict'], chart_json['samples'])
        if 'exp_venn' in chart_json:
            self.chart_exp_venn(chart_json['exp_venn'])
        if "exp_corr_file" in chart_json:
            group_dict = json.dumps(chart_json["group_dict"])
            exp_corr_file = chart_json["exp_corr_file"]
            exp_corr_tree_file = chart_json["exp_corr_tree_file"]
            self.chart_exp_corr(exp_corr_file, exp_corr_tree_file, group_dict)
        if "exp_pca_file" in chart_json:
            group_dict = chart_json["group_dict"]
            exp_pca_file = chart_json["exp_pca_file"]
            exp_pca_var_file = chart_json["exp_pca_var_file"]
            if "exp_pca_ellipse" in chart_json:
                exp_pca_ellipse = chart_json["exp_pca_ellipse"]
            else:
                exp_pca_ellipse = None
            self.chart_exp_pca(exp_pca_file, exp_pca_var_file, group_dict=group_dict, exp_pca_ellipse=exp_pca_ellipse, pcs=["PC1", "PC2"])
        if 'gene_diff_summary' and 'diff_method' in chart_json:
            self.chart_diffexp_stat(chart_json['gene_diff_summary'], cmp_list=chart_json['cmp_list'], soft=chart_json['diff_method'])
        if 'gene_diff_scatter' and 'diff_method' in chart_json:
            for each in chart_json["cmp_list"]:
                self.prok_diffexp_scatter(chart_json["gene_diff_scatter"].format(cmp1=each[0], cmp2=each[1]),
                                          '_vs_'.join(each), soft=chart_json['diff_method'])
        if 'srna_length' in chart_json:
            self.prok_srna_predict_length(chart_json['srna_length'])
        if 'srna_venn' in chart_json:
            self.prok_srna_stats_venn(chart_json['srna_venn'])
        if 'srna_rfam' in chart_json:
            self.prok_srna_rfam_pie(chart_json['srna_rfam'])
        if 'operon_xls' in chart_json:
            self.prok_operon_bar(chart_json['operon_xls'])
        if 'utr_xls' in chart_json:
            self.prok_utr_len_bar(chart_json['utr_xls'])
        if 'snp_anno' in chart_json:
            self.prok_snp_pie_bar(chart_json['snp_anno'], chart_json['samples'])
        if 'cluster_exp' and 'cluster_tree' and 'sample_tree' and 'subclusters' in chart_json:
            subcluster_list = glob.glob(chart_json['subclusters'])
            self.chart_geneset_cluster(chart_json['cluster_exp'], chart_json['cluster_tree'], chart_json['sample_tree'],
                                       subcluster_list, group_dict=chart_json["group_dict"],
                                       samples_order=chart_json['samples'])
        if 'geneset_detail' in chart_json:
            self.chart_geneset_venn_file(chart_json['geneset_detail'])
        if 'analysis_json' in chart_json:
            with open(chart_json['analysis_json'], 'r') as f:
                geneset_info = json.loads(f.read())
            geneset_results = geneset_info['geneset_results']
            if 'geneset_name' in geneset_info and 'diff_cog_class' in geneset_results:
                geneset_cog = os.path.join(geneset_results['diff_cog_class'], 'cog_class_table.xls')
                self.chart_geneset_class_cog(geneset_cog, [geneset_info['geneset_name']])
            if 'diff_go_class' in geneset_results:
                geneset_go_class = os.path.join(geneset_results['diff_go_class'], 'go_class_table.xls')
                self.prok_geneset_goclass(geneset_go_class)
            if 'diff_go_enrich' in geneset_results and 'geneset_name' in geneset_info:
                geneset_go_enrich = os.path.join(geneset_results['diff_go_enrich'], 'go_enrich_All_Diff_mRNA_gene.xls')
                self.chart_geneset_gorich(geneset_go_enrich, geneset_info['geneset_name'])
            if 'diff_kegg_class' in geneset_results and 'kegg_level' in chart_json:
                geneset_kegg_class = os.path.join(geneset_results['diff_kegg_class'], 'kegg_stat.xls')
                self.prok_geneset_kegg_bar(geneset_kegg_class, chart_json['kegg_level'])
            if 'geneset_name' and 'geneset_list' in geneset_info and 'diff_kegg_enrich' in geneset_results:
                kegg_enrich_filename = os.path.basename(geneset_info['geneset_list']) + '.DE.list.check.kegg_enrichment.xls'
                geneset_kegg_enrich = os.path.join(geneset_results['diff_kegg_enrich'], 'enrich', kegg_enrich_filename)
                self.chart_geneset_keggrich(geneset_kegg_enrich, geneset_info['geneset_name'])

        self.generate_html_sh()


if __name__ == '__main__':
    a = Chart()


    chart_json = {
        "samples": ["A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3"],
        # "group_dict": {"A":["A1","A2","A3"], "B": ["B1","B2","B3"], "C": ["C1","C2","C3"], "D": ["D1","D2","D3"]},
        "qc_file_before": "/mnt/lustre/users/sanger-dev/wpm2/workspace/20210825/Prokrna_smgo_vc6ehuord2lsghcif8n5v0/HiseqReadsStat/output/qualityStat/{sample_name}.l.qual_stat,/mnt/lustre/users/sanger-dev/wpm2/workspace/20210825/Prokrna_smgo_vc6ehuord2lsghcif8n5v0/HiseqReadsStat/output/qualityStat/{sample_name}.r.qual_stat",
        # "align_satu_r": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/saturation/satur_{sample_name}.eRPKM.xls.saturation.R",
        # "align_satu_p": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/saturation/satur_{sample_name}.eRPKM.xls.cluster_percent.xls",
        # "align_coverage": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/coverage/{sample_name}.geneBodyCoverage.txt",
        # "align_pos": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/distribution/{sample_name}.reads_distribution.txt",
        # "align_chr": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/MapAssessment/output/chr_stat/{sample_name}.bam_chr_stat.xls",
        # "assemble_step": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/RefrnaAssemble/RefassembleStat/output/trans_count_stat_200.txt",
        # "assemble_new": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/RefrnaAssemble/RefassembleStat/output/code_num.txt",
        # "gene_exp_ref": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/ref.gene.tpm.matrix",
        # "tran_exp_ref": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/ref.transcript.fpkm.matrix",
        # "gene_exp_all": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/gene.tpm.matrix",
        # "tran_exp_all": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/Quant/transcript.fpkm.matrix",
        # "annot_stat": "/mnt/ilustre/users/sanger-dev/workspace/20201116/Refrna_tsg_248805/AnnotMerge__1/output/allannot_class/all_stat.xls"
        }

    a.chart_json_batch(chart_json)
    # a.prok_chart_raw_qc(samples, qc_file_before)

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
