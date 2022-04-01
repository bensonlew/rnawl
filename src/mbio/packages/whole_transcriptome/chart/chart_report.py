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
                "110104": {
                    "main_coll": "qc_detail_long",
                    "name": "longRNA-seq碱基质量分布图",
                    "img": "*long_clean_qc_qual.box.pdf",
                    "category_name": "测序数据质控-longRNA-seq "
                },
                "110105": {
                    "main_coll": "qc_detail_long",
                    "name": "longRNA-seq碱基错误率分布图",
                    "img": "*long*_qc_error.line.pdf",
                    "category_name": "测序数据质控-longRNA-seq "
                },
                "110106": {
                    "main_coll": "qc_detail_long",
                    "name": "longRNA-seq碱基含量分布图",
                    "img": "*long*_qc_base.line.pdf",
                    "category_name": "测序数据质控-longRNA-seq "
                },
                "110201": {
                    "main_coll": "qc_detail_small",
                    "name": "smallRNA-seq碱基质量分布图",
                    "img": "*small_raw_qc_qual.box.pdf",
                    "category_name": "测序数据质控-smallRNA-seq "
                },
                "110202": {
                    "main_coll": "qc_detail_small",
                    "name": "smallRNA-seq碱基错误率分布图",
                    "img": "*small*_qc_error.line.pdf",
                    "category_name": "测序数据质控-smallRNA-seq "
                },
                "110601": {
                    "main_coll": "sg_transcripts_step",
                    "name": "转录本长度分布柱状图",
                    "img": "*assemble_length_distribution*.pdf",
                    "category_name": "转录本组装"
                },
                "110701": {
                    "main_coll": "sg_transcripts_seq_type",
                    "name": "新转录本类型分布饼图",
                    "img": "new_transcriptome.assemble_distribution*pdf",
                    "category_name": ""
                },
                "111704": {
                    "main_coll": "sg_specimen_after",
                    "name": "circRNA-seq碱基质量分布图",
                    "img": "*circle_clean_qc_qual.box.pdf",
                    "category_name": "测序数据质控-circRNA-seq "
                },
                "111705": {
                    "main_coll": "sg_specimen_after",
                    "name": "circRNA-seq碱基错误率分布图",
                    "img": "*circle_clean_qc_error.line.pdf",
                    "category_name": "测序数据质控-circRNA-seq "
                },
                "111706": {
                    "main_coll": "sg_specimen_after",
                    "name": "circRNA-seq碱基含量分布图",
                    "img": "*circle_clean_qc_base.line.pdf",
                    "category_name": "测序数据质控-circRNA-seq "
                },
                "130301": {
                    "main_coll": "sg_new_lncrna_predict_stat_detail",
                    "name": "新lncRNA预测Venn图",
                    "img": "*new.lncRNA_predict.venn.pdf",
                    "category_name": "lncRNA鉴定与预测"
                },
                "130401": {
                    "main_coll": "sg_lncrna_statistics_in_samples",
                    "name": "lncRNA统计图",
                    "img": "*sample_all*lnc_predict_stat*pdf",
                    "category_name": ""
                },
                "130501": {
                    "main_coll": "sg_lncrna_statistics_in_category",
                    "name": "lncRNA分类统计图",
                    "img": "*category_all*lnc_predict_stat*pdf",
                    "category_name": ""
                },
                "130901": {
                    "main_coll": "sg_mirna_stat_detail",
                    "name": "miRNA 统计图",
                    "img": "sample_all*mirna_predict_stat*pdf",
                    "category_name": "miRNA鉴定与预测"
                },
                "131001": {
                    "main_coll": "sg_srna_stat_detail",
                    "name": "sRNA 统计图",
                    "img": "*sRNAs_distribution.pie.pdf",
                    "category_name": ""
                },
                "131201": {
                    "main_coll": "sg_circrna_identify_in_category",
                    "name": "circRNA分类统计图",
                    "img": "all*circrna_predict_stat*pdf",
                    "category_name": "circRNA预测"
                },
                "131401": {
                    "main_coll": "sg_annotationstat_detail",
                    "name": "功能注释统计图-柱图",
                    "img": "*annot_*_stat*pdf",
                    "category_name": "转录组功能注释-功能注释统计"
                },
                "140201": {
                    "main_coll": "sg_expgraph_mrna",
                    "name": "表达量分布图-盒型图",
                    "img": "mRNA*box.pdf",
                    "category_name": "表达量分析"
                },
                "140301": {
                    "main_coll": "sg_expvenn_mrna",
                    "name": "样本间 Venn图",
                    "img": "*mRNA_*.venn.pdf",
                    "category_name": ""
                },
                "141301": {
                    "main_coll": "sg_expcorr_detail",
                    "name": "样本间相关性热图",
                    "img": "long_*.exp.heat_corr.pdf",
                    "category_name": ""
                },
                "141401": {
                    "main_coll": "sg_exppca_detail",
                    "name": "PCA图",
                    "img": "*exp_relation_pca*pdf",
                    "category_name": ""
                },
                "150201": {
                    "main_coll": "sg_diff_summary",
                    "name": "表达量差异统计柱状图",
                    "img": "*differential.summary*pdf",
                    "category_name": "表达量差异分析"
                },
                "150301": {
                    "main_coll": "sg_diff_graph_mrna",
                    "name": "表达量差异火山图",
                    "img": "*mRNA*diff.volcano*pdf",
                    "category_name": ""
                },
                "150302": {
                    "main_coll": "sg_diff_graph_mrna",
                    "name": "表达量差异散点图",
                    "img": "*mRNA*diff.scatter*pdf",
                    "category_name": ""
                },
                "170101": {
                    "main_coll": "sg_geneset_venn_class_detail",
                    "name": "Venn图",
                    "img": "diff_genesets.analysis*pdf",
                    "category_name": "基因集分析"
                },
                "170201": {
                    "main_coll": "sg_geneset_cluster_detail",
                    "name": "聚类热图",
                    "img": "*geneset.cluster.heat_corr.pdf",
                    "category_name": ""
                },
                "170301": {
                    "main_coll": "sg_geneset_cluster_detail_sub",
                    "name": "子聚类趋势图 ",
                    "img": "*cluster.*line.pdf",
                    "category_name": ""
                },
                "170401": {
                    "main_coll": "sg_geneset_cog_class_detail",
                    "name": "COG分类统计柱状图  ",
                    "img": "*cog_annot.gene_set.column.pdf",
                    "category_name": ""
                },
                "170501": {
                    "main_coll": "sg_geneset_go_class_detail",
                    "name": "GO分类统计柱形图",
                    "img": "*.go_annot.gene_set.column.pdf",
                    "category_name": ""
                },
                "170701": {
                    "main_coll": "sg_geneset_kegg_class_statistic",
                    "name": "Pathway分类统计柱状图",
                    "img": "*.kegg_annot.*.column.pdf",
                    "category_name": ""
                },
                "170802": {
                    "main_coll": "sg_geneset_go_enrich_detail",
                    "name": "GO富集分析结果图--柱形图(带折线）",
                    "img": "*.go_enrich.gene_set.*.pdf",
                    "category_name": ""
                },
                "170903": {
                    "main_coll": "sg_geneset_kegg_enrich_detail",
                    "name": "KEGG富集分析结果图--气泡图",
                    "img": "*kegg_enrich.gene_set.*.pdf",
                    "category_name": ""
                },
                "180201": {
                    "main_coll": "sg_snp_graphic",
                    "name": "不同区域分布饼图 ",
                    "img": "*snp.pos_stat.pie*pdf",
                    "category_name": "SNP/InDel 分析"
                },
                "180501": {
                    "main_coll": "sg_snp_graphic_type_stat",
                    "name": "SNP类型统计图--饼图",
                    "img": "*snp.type_stat.pie*pdf",
                    "category_name": ""
                },
                "180801": {
                    "main_coll": "sg_splicing_rmats_diff_stats",
                    "name": "差异可变剪切事件统计饼状图",
                    "img": "*diff_splice_stat*pie.pdf",
                    "category_name": ""
                },
                "181201": {
                    "main_coll": "sg_splicing_rmats_count_detail",
                    "name": "可变剪切事件统计图",
                    "img": "*all_*splice_stat*pdf",
                    "category_name": "可变剪切分析"
                },
                "181301": {
                    "main_coll": "sg_bias_first",
                    "name": "不同长度miRNA首位碱基偏好性分布图",
                    "img": "*first_bias_per*pdf",
                    "category_name": "miRNA结构分析"
                },
                "181401": {
                    "main_coll": "sg_bias_location",
                    "name": "miRNA不同位置碱基的偏好性分布图",
                    "img": "*all_loc_bias_per*pdf",
                    "category_name": "miRNA结构分析"
                },
                "181501": {
                    "main_coll": "sg_mirna_edit_detail",
                    "name": "miRNA碱基编辑类型分布图--堆积图",
                    "img": "*miRNAedit_distribution*pdf",
                    "category_name": ""
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
