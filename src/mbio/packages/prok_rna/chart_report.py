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
import cairosvg
from PIL import Image
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
    # geneset_name = os.path.basename(file_name).split(".")[0]
    geneset_name = 'All_Diff_mRNA'
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
            "110101": {
                "main_coll": "sg_specimen",
                "name": "质控数据统计表",
                "img": "*raw_qc_base.line.pdf",
                "category_name": "RawReads碱基组成分布图"
            },
            "110102": {
                "main_coll": "sg_specimen",
                "name": "质控数据统计表",
                "img": "*clean_qc_base.line.pdf",
                "category_name": "CleanReads碱基组成分布图"
            },
            "110103": {
                "main_coll": "sg_specimen",
                "name": "质控数据统计表",
                "img": "*raw_qc_error.line.pdf",
                "category_name": "RawReads碱基错误率分布图"
            },
            "110104": {
                "main_coll": "sg_specimen",
                "name": "质控数据统计表",
                "img": "*clean_qc_error.line.pdf",
                "category_name": "CleanReads碱基错误率分布图"
            },
            "110105": {
                "main_coll": "sg_specimen",
                "name": "质控数据统计表",
                "img": "*raw_qc_qual.box.pdf",
                "category_name": "RawReads碱基质量分布图"
            },
            "110106": {
                "main_coll": "sg_specimen",
                "name": "质控数据统计表",
                "img": "*clean_qc_qual.box.pdf",
                "category_name": "CleanReads碱基质量分布图"
            },
            "120201": {
                "main_coll": "sg_assessment_saturation",
                "name": "测序饱和度曲线表",
                "img": "*map_saturation.line.pdf",
                "category_name": "测序饱和度曲线图"
            },
            "120301": {
                "main_coll": "sg_assessment_coverage",
                "name": "基因覆盖度分布表",
                "img": "map_coverage.line.pdf",
                "category_name": "基因覆盖度分布图"
            },
            "130101": {
                "main_coll": "sg_annotation_stat",
                "name": "基础注释统计结果表",
                "img": "stats_annot_venn.venn.pdf",
                "category_name": "基础注释统计Venn图"
            },
            "130102": {
                "main_coll": "sg_annotation_stat",
                "name": "基础注释统计结果表",
                "img": "stats_annot_num_bar.bar.pdf",
                "category_name": "基础注释统计柱状图"
            },
            "131001": {
                "main_coll": "sg_annotation_cog",
                "name": "COG分类统计柱状表",
                "img": "annot_cog.cog_bar.pdf",
                "category_name": "COG分类统计柱状图"
            },
            "131501": {
                "main_coll": "sg_annotation_go",
                "name": "GO注释分类统计表",
                "img": "annot_go_pie.multi_pie.pdf",
                "category_name": "GO注释分类统计饼图"
            },
            "131502": {
                "main_coll": "sg_annotation_go",
                "name": "GO注释分类统计表",
                "img": "annot_go_bar.multi_bar.pdf",
                "category_name": "GO注释分类统计柱状图"
            },
            "131601": {
                "main_coll": "sg_annotation_kegg",
                "name": "Pathway信号通路注释结果统计表",
                "img": "annot_kegg_bar.kegg_bar.pdf",
                "category_name": "Pathway分类统计柱状图"
            },
            "140201": {
                "main_coll": "sg_exp_graph",
                "name": "表达量图形表",
                "img": "*exp_distribution.box.pdf",
                "category_name": "表达量分布盒形图"
            },
            "140202": {
                "main_coll": "sg_exp_graph",
                "name": "表达量图形表",
                "img": "*exp_distribution.violin.pdf",
                "category_name": "表达量分布小提琴图"
            },
            "140203": {
                "main_coll": "sg_exp_graph",
                "name": "表达量图形表",
                "img": "*exp_distribution.density.pdf",
                "category_name": "表达量分布密度图"
            },
            "140301": {
                "main_coll": "sg_exp_venn",
                "name": "表达量venn表",
                "img": "all.exp.venn.pdf",
                "category_name": "表达量venn图"
            },
            "140401": {
                "main_coll": "sg_exp_corr",
                "name": "样本间相关性系数表",
                "img": "exp.heatmap.heat_corr.pdf",
                "category_name": "样本间相关性热图"
            },
            "140501": {
                "main_coll": "sg_exp_pca",
                "name": "主成分解释表",
                "img": "all.exp_relation_pca.scatter.pdf",
                "category_name": "样本间PCA图"
            },
            "150201": {
                "main_coll": "sg_diff",
                "name": "表达量差异统计表",
                "img": "*.diffexp_summary.bar.pdf",
                "category_name": "表达量差异统计图-柱状图"
            },
            "150202": {
                "main_coll": "sg_diff",
                "name": "表达量差异统计表",
                "img": "*.diffexp_summary.stacked_bar.pdf",
                "category_name": "表达量差异统计图-堆积图"
            },
            "150301": {
                "main_coll": "sg_diff",
                "name": "表达量差异火山表",
                "img": "*.diffexp.volcano.pdf",
                "category_name": "表达量差异火山图"
            },
            "150401": {
                "main_coll": "sg_diff",
                "name": "表达量差异散点表",
                "img": "*.diffexp.scatter.pdf",
                "category_name": "表达量差异散点图"
            },
            "160101": {
                "main_coll": "sg_geneset_venn",
                "name": "基因集venn表",
                "img": "geneset.venn.pdf",
                "category_name": "基因集venn图"
            },
            "160201": {
                "main_coll": "sg_geneset_cluster",
                "name": "heatmap分析表",
                "img": "geneset.cluster.heat_corr.pdf",
                "category_name": "基因集聚类分析热图"
            },
            "160301": {
                "main_coll": "sg_geneset_cluster",
                "name": "子聚类heatmap分析表",
                "img": "subcluster*pdf",
                "category_name": "基因集子聚类趋势图"
            },
            "160401": {
                "main_coll": "sg_geneset_cog_class",
                "name": "基因集COG分类统计表",
                "img": "cog_annot.gene_set.column.pdf",
                "category_name": "基因集COG分类统计柱状图"
            },
            "160501": {
                "main_coll": "sg_geneset_go_class",
                "name": "基因集GO分类统计表",
                "img": "*.go_bar_geneset.go_bar.pdf",
                "category_name": "基因集GO注释柱形图"
            },
            "160601": {
                "main_coll": "sg_geneset_kegg_class",
                "name": "基因集KeggPathway统计表",
                "img": "*annot_kegg_bar.kegg_bar.pdf",
                "category_name": "基因集KeggPathway分类统计柱状图"
            },
            "160701": {
                "main_coll": "sg_geneset_go_enrich",
                "name": "GO富集分析统计表",
                "img": "*_go_enrich_bar.shadowbar.pdf",
                "category_name": "GO富集分析|柱形图",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "160702": {
                "main_coll": "sg_geneset_go_enrich",
                "name": "GO富集分析统计表",
                "img": "*_all_go_enrich_dense_bubble.densebubble.pdf",
                "category_name": "GO富集分析|气泡图",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "160703": {
                "main_coll": "sg_geneset_go_enrich",
                "name": "GO富集分析统计表",
                "img": "*_all_go_enrich_scatter_bubble.scatterbubble.pdf",
                "category_name": "GO富集分析|分散型气泡图",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "160801": {
                "main_coll": "sg_geneset_kegg_enrich",
                "name": "KEGG富集分析统计表",
                "img": "*_enrichkegg.shadowbar.pdf",
                "category_name": "KEGG富集分析|柱形图",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "160802": {
                "main_coll": "sg_geneset_kegg_enrich",
                "name": "KEGG富集分析统计表",
                "img": "*_enrichkegg.densebubble.pdf",
                "category_name": "KEGG富集分析|气泡图",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "162001": {
                "main_coll": "sg_geneset_ipath",
                "name": "代谢通路数据",
                "img": "../DiffGenesetAllPipline/Ipath/Metabolic_pathways.svg",
                "category_name": "iPath代谢通路图"
            },
            "170201": {
                "main_coll": "sg_srna_predict",
                "name": "sRNA长度统计表",
                "img": "srna_length_bar.bar.pdf",
                "category_name": "sRNA长度统计"
            },
            "170401": {
                "main_coll": "sg_srna_anno_venn",
                "name": "注释结果Venn统计表",
                "img": "srna_annot_venn.venn.pdf",
                "category_name": "注释结果Venn统计"
            },
            "170601": {
                "main_coll": "sg_srna_anno",
                "name": "Rfam注释统计表",
                "img": "srna_rfam_pie.pie.pdf",
                "category_name": "Rfam注释统计图"
            },
            "180101": {
                "main_coll": "sg_structure_operon",
                "name": "操纵子长度统计表",
                "img": "operon_len_bar.bar.pdf",
                "category_name": "操纵子长度统计图"
            },
            "180201": {
                "main_coll": "sg_structure_operon",
                "name": "操纵子基因数量统计表",
                "img": "operon_gene_num_bar.bar.pdf",
                "category_name": "操纵子基因数量统计图"
            },
            "180601": {
                "main_coll": "sg_structure_utr",
                "name": "UTR长度分布表",
                "img": "UTR*_len_bar.bar.pdf",
                "category_name": "UTR长度分布图"
            },
            "181001": {
                "main_coll": "sg_snp",
                "name": "SNPIndel不同区域统计表",
                "img": "*snp_regions_pie.pie.pdf",
                "category_name": "区域分布饼图"
            },
            "181101": {
                "main_coll": "sg_snp",
                "name": "SNP深度统计表",
                "img": "*depth_stats_pie.pie.pdf",
                "category_name": "SNP深度统计饼图"
            },
            "181102": {
                "main_coll": "sg_snp",
                "name": "SNP深度统计表",
                "img": "*depth_stats_bar.bar.pdf",
                "category_name": "SNP深度统计柱状图"
            },
            "181201": {
                "main_coll": "sg_snp",
                "name": "SNP频率统计表",
                "img": "*frequency_stats_pie.pie.pdf",
                "category_name": "SNP频率统计饼图"
            },
            "181202": {
                "main_coll": "sg_snp",
                "name": "SNP频率统计表",
                "img": "*frequency_stats_bar.bar.pdf",
                "category_name": "SNP频率统计柱状图"
            },
            "181301": {
                "main_coll": "sg_snp",
                "name": "SNP类型统计表",
                "img": "*type_stats_pie.pie.pdf",
                "category_name": "SNP类型统计饼图"
            },
            "181302": {
                "main_coll": "sg_snp",
                "name": "SNP类型统计表",
                "img": "*type_stats_bar.bar.pdf",
                "category_name": "SNP类型统计柱状图"
            },
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
        if str(pdf_choose).endswith('pdf'):
            png_path = os.path.splitext(os.path.basename(pdf_choose))[0]
            # command = "gs -o png/{}.jpg -sDEVICE=jpeg  -r256 {} && convert png/{}.jpg png/{}.png".format(png_path, pdf_choose, png_path, png_path)
            pdf2png = "gs -o png/{}.jpg -sDEVICE=jpeg -r256 {}".format(png_path, pdf_choose)
            os.system(pdf2png)

            # work out dynamic dimensions of png
            img = Image.open('png/{}.jpg'.format(png_path))
            width, height = img.size
            resize = '{}x{}'.format(str(int(width/3.6)), str(int(height/3.6)))

            jpg2png = "convert -resize {} png/{}.jpg png/{}.png".format(resize, png_path, png_path)
            os.system(jpg2png)
        elif str(pdf_choose).endswith('svg'):
            png_path = os.path.splitext(os.path.basename(pdf_choose))[0]
            cairosvg.svg2png(url=pdf_choose, write_to=os.path.join("./png/", png_path + '.png'))
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
