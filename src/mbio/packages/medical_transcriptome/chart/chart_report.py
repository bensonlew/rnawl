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
reload(sys)
sys.setdefaultencoding('utf-8')


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
            "1001": {
                "categroy_name": "表达量差异统计",
                "img": u"../DiffexpBatch/all.diffexp_summary.bar2.pdf",
                "name": u"表达量统计图-堆积图",
                "main_coll": u"sg_diff"
            },
            "1002": {
                "categroy_name": u"表达量差异统计",
                "img": u"../DiffexpBatch/*.diffexp.volcano.pdf",
                "name": u"表达量差异火山图",
                "main_coll": u"sg_diff"
            },
            "1003": {
                "categroy_name": u"表达量差异统计",
                "img": u"../DiffexpBatch/*.diffexp_ma.scatter.pdf",
                "name": u"表达量差异MA图",
                "main_coll": u"sg_diff"
            },
            "1004": {
                "categroy_name": u"差异基因Venn分析",
                "img": u"diff_genesets.analysis.venn.pdf",
                "name": u"Venn图",
                "main_coll": u"sg_diff_geneset_venn"
            },
            "1005": {
                "categroy_name": u"差异基因聚类分析",
                "img": u"geneset.cluster.heat_corr.pdf",
                "name": u"聚类热图",
                "main_coll": u"sg_diff_geneset_cluster"
            },
            "1007": {
                "categroy_name": u"差异基因功能注释分析",
                "img": u"*.go_annot.gene_set.column.pdf",
                "name": u"GO分类统计图-柱状图",
                "main_coll": u"sg_diff_geneset_go_class",
                "main_dict": {"gene_set": u"get_geneset_by_file_name"}
            },
            "1008": {
                "categroy_name": u"差异基因功能注释分析",
                "img": u"*.kegg_annot.gene_set.column.pdf",
                "name": u"KEGG分类统计图-柱状图",
                "main_coll": u"sg_diff_geneset_kegg_class",
                "main_dict": {"gene_set": u"get_geneset_by_file_name"}
            },
            "1009": {
                "categroy_name": u"差异基因功能注释分析",
                "img": u"*.reactome_annot.gene_set.column.pdf",
                "name": u"Reactome分类统计图-柱状图",
                "main_coll": u"sg_diff_geneset_reactome_class",
                "main_dict": {"gene_set": u"get_geneset_by_file_name"}
            },
            "1010": {
                "categroy_name": u"差异基因功能注释分析",
                "img": "*.do_annot.gene_set.column.pdf",
                "name": u"DO分类统计图-柱状图",
                "main_coll": "sg_diff_geneset_do_class",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "1011": {
                "categroy_name": u"差异基因功能富集分析",
                "img": "*.go_enrich.gene_set.bar_line.pdf",
                "name": u"GO富集分析结果图-柱形图（带折线）",
                "main_coll": "sg_diff_geneset_go_enrich",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "1012": {
                "categroy_name": u"差异基因功能富集分析",
                "img": "*.kegg_enrich.gene_set.buble.pdf",
                "name": u"KEGG富集分析结果图-气泡图",
                "choose_func": "choose_js_large",
                "main_coll": "sg_diff_geneset_kegg_enrich",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "1013": {
                "categroy_name": u"差异基因功能富集分析",
                "img": "*.reactome_enrich.gene_set.buble.pdf",
                "name": u"Reactome富集分析结果图-气泡图",
                "main_coll": "sg_diff_geneset_reactome_enrich",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "1014": {
                "categroy_name": u"差异基因功能富集分析",
                "img": "*.do_enrich.gene_set.bar_line.pdf",
                "name": u"DO富集分析结果图-柱形图（带折线）",
                "main_coll": "sg_diff_geneset_do_enrich",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}

            },
            "1015": {
                "categroy_name": u"样本关系分析",
                "img": "all.exp.venn.pdf",
                "name": u"样本间 Venn图",
                "main_coll": "sg_geneset_venn"
            },
            "1016": {
                "categroy_name": u"样本关系分析",
                "img": "all.exp.heat_corr.pdf",
                "name": u"样本间相关性热图",
                "main_coll": "sg_exp_corr"
            },
            "1017": {
                "categroy_name": u"样本关系分析",
                "img": "all.exp_relation_pca.scatter.pdf",
                "name": u"PCA图",
                "main_coll": "sg_exp_pca"
            },
            "1018": {
                "categroy_name": u"差异可变剪切分析",
                "img": "*.diff_splice_stat.pie.pdf",
                "name": u"差异可变剪切事件统计图",
                "main_coll": "sg_splicing_rmats"
            },
            "101802": {
                "categroy_name": u"差异可变剪切分析",
                "img": "*.diff_splice_stat.column.pdf",
                "name": u"差异可变剪切模式变化统计图",
                "main_coll": "sg_splicing_rmats"
            },
            "1019": {
                "categroy_name": u"基因融合分析",
                "img": "*gene_fusion.pdf",
                "name": u"基因融合circos图",
                "main_coll": "sg_gene_fusion"

            },
            "1020": {
                "categroy_name": u"SNP/InDel分析",
                "img": "*.snp.pos_stat.pie.pdf",
                "name": u"不同区域分布饼图",
                "main_coll": "sg_snp",
                "main_dict": {"sample": "get_sample_by_file_name"}
            },
            "102002": {
                "categroy_name": u"SNP/InDel分析",
                "img": "*.snp.type_stat.column.pdf",
                "name": u"SNP类型统计柱状图",
                "main_coll": "sg_snp",
                "main_dict": {"sample": "get_sample_by_file_name"}
            },
            "1021": {
                "categroy_name": u"测序数据质控-质控数据统计",
                "img": "*.clean_qc_error.line.pdf",
                "name": u"碱基错误率分布图",
                "main_coll": "sg_qc",
                "main_dict": {"sample": "get_sample_by_file_name"}
            },
            "1022": {
                "categroy_name": u"测序数据质控-质控数据统计",
                "img": "*.raw_qc_base.line.pdf",
                "name": u"碱基含量分布图",
                "main_coll": "sg_qc",
                "main_dict": {"sample": "get_sample_by_file_name"}
            },
            "1023": {
                "categroy_name": u"测序数据质控-质控数据统计",
                "img": "*.raw_qc_qual.box.pdf",
                "name": u"盒形图",
                "main_coll": "sg_qc",
                "main_dict": {"sample": "get_sample_by_file_name"}
            },
            "1024": {
                "categroy_name": u"序列比对分析-转录组质量评估",
                "img": "*.align_coverage.line.pdf",
                "name": u"测序覆盖度分布图",
                "main_coll": "sg_assessment_coverage",
                "main_dict": {"sample": "get_sample_by_file_name"}
            },
            "1025": {
                "categroy_name": u"序列比对分析-转录组质量评估",
                "img": "*.align_pos_dist.pie.pdf",
                "name": u"不同区域Reads分布统计饼图",
                "main_coll": "sg_assessment_distribution",
                "main_dict": {"sample": "get_sample_by_file_name"}
            },
            "1026": {
                "categroy_name": u"表达量分布",
                "img": "G_samples.exp_distribution.box.pdf",
                "name": u"表达量分布-盒型图",
                "main_coll": "sg_exp_graph"
            }
        }

    def run(self, config="report_config.json"):
        self.filter_chart()
        # self.pdf2png()
        self.write_config(config)

    def write_config(self, config):
        with open(config, 'w') as fw:
            fw.write(json.dumps(self.report_config_out, indent=4, ensure_ascii=False))

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
                "categroy_name": img_dict["categroy_name"],
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
