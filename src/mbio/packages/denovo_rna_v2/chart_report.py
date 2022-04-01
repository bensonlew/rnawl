# coding=utf-8
from __future__ import division
import xml.etree.ElementTree as ET
from biocluster.config import Config
import re
import collections
import json
import subprocess
import gridfs
import glob
import os
import sys
import pandas as pd
import json
import math
import numpy as np
from collections import OrderedDict


def choose_js_large(files):
    '''
    选择["dataset"][0]["source"][0]中最多的文件
    '''
    file_list = list()
    for file_ in files:
        file_pre = os.path.splitext(file_)
        file_js = file_pre + ".js"
        json_d = get_json_from_js(file_js)
        file_list.append([file_pre, file_js, len(json_d["dataset"][0]["source"][0])])

    file_list_sort = sorted(file_list, key=lambda x:x[2], reverse=True)
    return file_list[0][0]

def get_json_from_js(js_file):
    with open(js_file, 'r') as f:
        b = f.read()
        json_d = json.loads(b[14:])
    return json_d

def get_geneset_by_file_name(file_name):
    geneset_name = os.path.basename(file_name).split(".")[0]
    return geneset_name

def get_sample_by_file_name(file_name):
    sample_name = os.path.basename(file_name).split(".")[0]
    return sample_name


class ChartReport(object):
    def __init__(self):
        """
        设置数据库，连接到mongod数据库，kegg_ko,kegg_gene,kegg_pathway_png三个collections2
        """

        """
        class: sample[A1, B1], group[A, B], classify: [samples, groups], level[G, T], kind[ref, new, all], category[mRNA, lncRNA]
        """
        super(ChartReport, self).__init__()
        self.report_init()
        os.system("mkdir -p png")



    def report_init(self):
        '''
        "categroy_name": 目录
        "img": 静态图片名称, 通过glob.glob 匹配图片
        "name": "表达量统计图-堆积图",
        "main_coll": 主表collection名称
        "choose_func": 图片存在多张时通过该函数筛选图片
        "choose_main"
        '''
        self.convert_list = list()
        self.report_config_out = dict()
        self.report_config = {
            "110301": {
                "main_coll": "sg_specimen_after",
                "name": "碱基质量分布图",
                "img": "*.clean_qc_qual.box.pdf",
                "category_name": "测序数据质控-质控数据统计 "
            },
            "110302": {
                "main_coll": "sg_specimen_after",
                "name": "碱基错误率分布图",
                "img": "*.clean_qc_error.line.pdf",
                "category_name": "测序数据质控-质控数据统计 "
            },
            "110303": {
                "main_coll": "sg_specimen_after",
                "name": "碱基含量分布图",
                "img": "*.clean_qc_base.line.pdf",
                "category_name": "测序数据质控-质控数据统计 "
            },
            "120301": {
                "main_coll": "sg_assembly_len",
                "name": "转录本长度分布图",
                "img": "Transcript.length_distribution.columns.pdf",
                "category_name": "转录组从头组装-转录本长度分布"
            },
            "130101": {
                "main_coll": "sg_annotationstat",
                "name": "功能注释统计图-柱图",
                "img": "annot_gene_stat.column.pdf",
                "category_name": "转录组功能注释-功能注释统计"
            },
            "130301": {
                "main_coll": "sg_annotationnr",
                "name": "NR注释物种分布饼图",
                "img": "Gene.nr_species.pie.pdf",
                "category_name": "转录组功能注释-NR"
            },
            "130901": {
                "main_coll": "sg_annotationpfam",
                "name": "Pfam注释柱状图",
                "img": "Pfam.annotation_of_Gene.columns.pdf",
                "category_name": "转录组功能注释-Pfam"
            },
            "131001": {
                "main_coll": "sg_annotationcog",
                "name": "COG分类统计柱状图",
                "img": "COG_Gene.function_classification.columns.pdf",
                "category_name": "转录组功能注释-COG"
            },
            "131101": {
                "main_coll": "sg_annotationgo",
                "name": "GO分类统计图-BP",
                "img": "Gene_GO_level2.classification_BP.pie.pdf",
                "category_name": "转录组功能注释-GO"
            },
            "131201": {
                "main_coll": "sg_annotationkegg",
                "name": "Pathway分类统计柱状图",
                "img": "Gene_Histogram_of_KEGG.barline.barline.pdf",
                "category_name": "转录组功能注释-KEGG"
            },
            "150201": {
                "main_coll": "sg_tf",
                "name": "转录因子家族统计图",
                "img": "Unigene_TF.family.columns.pdf",
                "category_name": "转录因子家族分析-转录因子家族分析"
            },
            "160201": {
                "main_coll": "sg_expgraph",
                "name": "表达量分布图-盒型图",
                "img": "Genes_samples_exp.distribution.box.box.pdf",
                "category_name": "表达量统计- 表达量分布"
            },
            "160301": {
                "main_coll": "sg_expvenn",
                "name": "样本间 Venn图",
                "img": "all.exp.venn.pdf",
                "category_name": "样本关系分析"
            },

            "160401": {
                "main_coll": "sg_expcorr",
                "name": "样本间相关性热图",
                "img": "all.exp.heat_corr.pdf",
                "category_name": "样本关系分析"
            },
            "160501": {
                "main_coll": "sg_exppca",
                "name": "PCA图",
                "img": "all.exp_relation_pca.scatter.pdf",
                "category_name": "样本关系分析"
            },
            "170301": {
                "main_coll": "sg_diff",
                "name": "表达量差异统计柱状图",
                "img": "Gene.differential.summary.bar_h.bar.pdf",
                "category_name": "表达量差异统计"
            },
            "170401": {
                "main_coll": "diff_graph",
                "name": "表达量差异火山图",
                "img": "*.diff.volcano.volcano.pdf",
                "category_name": "表达量差异统计"
            },
            "170402": {
                "main_coll": "diff_graph",
                "name": "表达量差异散点图",
                "img": "*.diff.scatter.scatter.pdf",
                "category_name": "表达量差异统计"
            },
            "180101": {
                "main_coll": "sg_geneset_venn",
                "name": "Venn图",
                "img": "diff_genesets.analysis.venn.pdf",
                "category_name": "基因集分析-Venn分析"
            },
            "180201": {
                "main_coll": "sg_geneset_cluster",
                "name": "聚类热图",
                "img": "geneset.cluster.heat_corr.pdf",
                "category_name": "基因集分析-聚类分析"
            },
            "180301": {
                "main_coll": "sg_geneset_cluster",
                "name": "子聚类趋势图 ",
                "img": "cluster.*.line.pdf",
                "category_name": "基因集分析-聚类分析"
            },
            "180401": {
                "main_coll": "sg_geneset_cog_class",
                "name": "COG分类统计柱状图  ",
                "img": "*.cog_annot.gene_set.column.pdf",
                "category_name": "基因集功能注释--COG",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "180501": {
                "main_coll": "sg_geneset_go_class",
                "name": "GO分类统计柱形图",
                "img": "*.go_annot.gene_set.column.pdf",
                "category_name": "基因集功能注释--GO",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "180701": {
                "main_coll": "sg_geneset_go_enrich",
                "name": "GO富集分析结果图--柱形图(带折线）",
                "img": "*.go_enrich.gene_set.bar_line.pdf",
                "category_name": "基因集功能富集--GO",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "180803": {
                "main_coll": "sg_geneset_kegg_enrich",
                "name": "KEGG富集分析结果图--气泡图",
                "img": "*.kegg_enrich.gene_setgo_enrich_buble1.pdf",
                "category_name": "基因集功能富集--KEGG",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "180901": {
                "main_coll": "sg_geneset_kegg_class",
                "name": "Pathway分类统计柱状图",
                "img": "*.kegg_annot.*.column.pdf",
                "category_name": "基因集功能注释--KEGG",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "190101": {
                "main_coll": "sg_cds",
                "name": "CDS长度分布图",
                "img": "Gene.length_distribution_of_CDS.columns.pdf",
                "category_name": "基因结构分析--CDS预测"
            },
            "190401": {
                "main_coll": "sg_snp",
                "name": "SNP位点统计图-柱状图",
                "img": "*.snp.pos_stat.pie.pdf",
                "category_name": "基因结构分析--SNP 分析"
            },
            "190502": {
                "main_coll": "sg_snp",
                "name": "SNP测序深度统计图-饼图",
                "img": "*.snp.pos_stat.pie.pdf",
                "category_name": "基因结构分析--SNP 分析"
            },
            "190701": {
                "main_coll": "sg_ssr",
                "name": "SSR统计图",
                "img": "SSR.stat.pie.pdf",
                "category_name": "基因结构分析--SSR分析"
            },
            "190702": {
                "main_coll": "sg_ssr",
                "name": "SSR类型分布统计图",
                "img": "",
                "category_name": "基因结构分析--SSR分析"
            }
        }

    def run(self, config="report_config.json"):
        self.filter_chart()
        # self.pdf2png()
        self.write_config(config)

    def write_config(self, config):
        with open(config, 'w') as fw:
            fw.write(json.dumps(self.report_config_out, indent=4, ensure_ascii=False).encode('utf8').decode())

    def filter_chart(self):
        '''
        过滤图片
        '''
        for img_id, img_dict in self.report_config.items():
            files = glob.glob(img_dict["img"])
            if len(files) > 1:
                pdf_choose = self.choose(files, img_dict)
            elif len(files) == 1:
                pdf_choose = files[0]
            else:
                print("该分析未找到结果图片 {}".format(img_dict))
                continue

            png_path = self.pdf2png(pdf_choose)

            self.report_config_out[img_id] = {
                "img": png_path + ".png",
                "pdf_choose": pdf_choose,
                "name": img_dict["name"],
                "category_name": img_dict["category_name"],
                "main_coll": img_dict["main_coll"]
            }

            if "main_dict" in img_dict:
                main_dict_out = self.convert_main_dict(img_dict["main_dict"], pdf_choose)
                self.report_config_out[img_id].update({
                    "main_dict": main_dict_out
                })


    def convert_main_dict(self, main_dict, pdf_choose):
        '''
        获取主表筛选字段
        '''
        main_dict_choosed = dict()
        for k, v in main_dict.items():
            func = globals()[v]
            main_dict_choosed[k] = func(pdf_choose)
        return main_dict_choosed


    def choose(self, files, img_dict):
        '''
        根据配置文件函数筛选图片
        '''
        if "choose_func" in img_dict:
            func_name = img_dict["choose_func"]
            func = globals()[func_name]
            choosed = files[0]
            try:
                choosed = func(files)
            except Exception as e:
                print("筛选图片错误")
                print(e)
            return choosed
        else:
            return files[0]

    def pdf2png(self, pdf_choose):
        png_path = os.path.splitext(os.path.basename(pdf_choose))[0]
        command = "gs -o png/{}.jpg -sDEVICE=jpeg  -r256 {} && convert png/{}.jpg png/{}.png".format(png_path, pdf_choose, png_path, png_path)
        os.system(command)
        return png_path

    # def pdf2png(self):
    #     os.system("mkdir -p png")
    #     for img_id, pdf_choose in self.convert_list:
    #         png_path = os.path.splitext(os.path.basename(pdf_choose))[0]
    #         command = "gs -o png/{}.jpg -sDEVICE=jpeg  -r256 {} && convert png/{}.jpg png/{}.png".format(png_path, pdf_choose, png_path, png_path)
    #         print(command)
    #         os.system(command)
    #         self.report_config_out[img_id] = {
    #             "img": png_path + ".png",
    #             ""
    #         }



if __name__ == '__main__':
    a = ChartReport()
    a.run()
