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
            "120201": {
                "main_coll": "sg_assessment_distribution",
                "name": "不同区域Reads分布统计饼图 ",
                "img": "*.align_pos_dist.pie.pdf",
                "category_name": "序列比对分析-转录组质量评估"
            },
            "120301": {
                "main_coll": "sg_assessment_chrom_distribution",
                "name": "不同染色体Reads分布统计柱状图 ",
                "img": "*.align_chr_dist.bar.pdf",
                "category_name": "序列比对分析-转录组质量评估"
            },
            "120401": {
                "main_coll": "sg_assessment_saturation",
                "name": "测序饱和度曲线图",
                "img": "*.align.satu.line.pdf",
                "category_name": "序列比对分析-转录组质量评估"
            },
            "120501": {
                "main_coll": "sg_assessment_coverage",
                "name": "测序覆盖度分布图",
                "img": "all.align_coverage.line.pdf",
                "category_name": "序列比对分析-转录组质量评估"
            },
            "140102": {
                "main_coll": "sg_annotationstat",
                "name": "功能注释统计图-柱图 ",
                "img": "annot_gene_stat.column.pdf",
                "category_name": "功能注释与查询-功能注释统计"
            },
            "150201": {
                "main_coll": "sg_expgraph",
                "name": "表达量分布图-盒型图",
                "img": "G_samples.exp_distribution.box.pdf",
                "category_name": "表达量统计- 表达量分布"
            },
            "150301": {
                "main_coll": "sg_expvenn",
                "name": "样本间 Venn图",
                "img": "all.exp.venn.pdf",
                "category_name": "样本关系分析"
            },
            "150401": {
                "main_coll": "sg_expcorr",
                "name": "样本间相关性热图",
                "img": "all.exp.heat_corr.pdf",
                "category_name": "样本关系分析"
            },
            "150501": {
                "main_coll": "sg_expcorr",
                "name": "样本间pca图",
                "img": "all.exp_relation_pca.scatter.pdf",
                "category_name": "样本关系分析"
            },
            "160301": {
                "main_coll": "sg_diff",
                "name": "表达量差异火山图",
                "img": "../DiffexpBatch/*.diffexp.volcano.pdf",
                "category_name": "表达量差异统计"
            },
            "160401": {
                "main_coll": "sg_diff",
                "name": "表达量差异散点图",
                "img": "../DiffexpBatch/*.diffexp.scatter.pdf",
                "category_name": "表达量差异统计"
            },
            "170101": {
                "main_coll": "sg_geneset_venn",
                "name": "Venn图",
                "img": "diff_genesets.analysis.venn.pdf",
                "category_name": "Venn分析"
            },
            "170201": {
                "main_coll": "sg_geneset_cluster",
                "name": "聚类热图",
                "img": "geneset.cluster.heat_corr.pdf",
                "category_name": "聚类分析"
            },
            "170401": {
                "main_coll": "sg_geneset_cog_class",
                "name": "COG分类统计柱状图  ",
                "img": "*cog_annot.gene_set.column.pdf",
                "category_name": "功能注释--COG"
            },
            "170501": {
                "main_coll": "sg_geneset_go_class",
                "name": "GO分类统计柱形图",
                "img": "*go_annot.gene_set.column.pdf",
                "category_name": "功能注释--GO",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "170704": {
                "main_coll": "sg_geneset_go_enrich",
                "name": "GO富集分析结果图--柱形图(带折线）",
                "img": "*.go_enrich.gene_set.bar_line.pdf",
                "category_name": "功能富集--GO",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "170802": {
                "main_coll": "sg_geneset_kegg_enrich",
                "name": "KEGG富集分析结果图--气泡图",
                "img": "*.kegg_enrich.gene_set.buble.pdf",
                "category_name": "功能富集--KEGG",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "170901": {
                "main_coll": "sg_geneset_kegg_class",
                "name": "Pathway分类统计柱状图",
                "img": "*.kegg_annot.*.column.pdf",
                "category_name": "功能注释--KEGG",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "180201": {
                "main_coll": "sg_snp",
                "name": "不同区域分布饼图 ",
                "img": "*.snp.pos_stat.pie.pdf",
                "category_name": "SNP/Indel 分析-SNP/InDel不同区域分布"
            },
            "180501": {
                "main_coll": "sg_snp",
                "name": "SNP类型统计图--饼图 ",
                "img": "*.snp.type_stat.pie.pdf",
                "category_name": "SNP/Indel 分析-SNP统计"
            },
            "181201": {
                "main_coll": "sg_splicing_rmats_count",
                "name": "可变剪切事件统计图",
                "img": "all_JCsplice_stat.column.pdf",
                "category_name": "可变剪切分析-可变剪切事件统计"
            },
            "181202": {
                "main_coll": "sg_splicing_rmats_count",
                "name": "单样本可变剪切事件统计图",
                "img": "",
                "category_name": "可变剪切分析-可变剪切事件统计"
            },
            "210402": {
                "main_coll": "sg_exp_batch",
                "name": "PCA图",
                "img": "all.exp_relation_pca.scatter.pdf",
                "category_name": "样本关系分析"
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
            print img_dict
            self.report_config_out[img_id] = {
                "img": png_path + ".png",
                "pdf_choose": pdf_choose,
                "name": img_dict["name"],
                "categroy_name": img_dict["category_name"],
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
