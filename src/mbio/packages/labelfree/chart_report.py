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
import cairosvg


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
    geneset_name = os.path.basename(file_name).split("____")[0]
    if "up" not in geneset_name:
        if "down" not in geneset_name:
            if "all" not in geneset_name:
                geneset_name = geneset_name + "_all"
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
            "110201": {
                "main_coll": "sg_protein_info",
                "name": "蛋白质分子量分布图",
                "img": "./output/mw.bar.pdf",
                "category_name": "数据与质控评估"
            },
            "110202": {
                "main_coll": "sg_protein_info",
                "name": "蛋白质覆盖度分布图",
                "img": "./output/seq_cov.bar.pdf",
                "category_name": "数据与质控评估"
            },
            "110203": {
                "main_coll": "sg_protein_info",
                "name": "蛋白质信息统计图",
                "img": "./output/info.bar.pdf",
                "category_name": "数据与质控评估"
            },
            "110301": {
                "main_coll": "sg_peptide_num",
                "name": "肽段数量分布图",
                "img": "./output/pep_num.bar.pdf",
                "category_name": "数据与质控评估"
            },
            "110302": {
                "main_coll": "sg_peptide_error",
                "name": "肽段误差分布图",
                "img": "../DrawStaticplot/output/dMass.pdf",
                "category_name": "数据与质控评估"
            },
            "110303": {
                "main_coll": "sg_peptide_len",
                "name": "肽段长度分布",
                "img": "./output/pep_len.bar.pdf",
                "category_name": "数据与质控评估"
            },
            "120101": {
                "main_coll": "sg_express_corr",
                "name": "相关性热图",
                "img": "./output/sam_corr.heatmap.pdf",
                "category_name": "样本相关性分析"
            },
            "120201": {
                "main_coll": "sg_express_pca",
                "name": "PCA分析图",
                "img": "./output/sam_pca.scatter.pdf",
                "category_name": "样本相关性分析"
            },
            "130101": {
                "main_coll": "sg_annotation_stat",
                "name": "注释结果柱状图",
                "img": "./output/anno_stat.bar.pdf",
                "category_name": "全蛋白功能注释"
            },
            "130201": {
                "main_coll": "sg_annotation_go",
                "name": "GO分类统计直方图",
                "img": "./output/go_lev2_bar.go_bar.pdf",
                "category_name": "全蛋白功能注释"
            },
            # "130202": {
            #     "main_coll": "sg_annotation_go",
            #     "name": "GO分类统计饼图",
            #     "img": "./output/go_lev2_pie.go_pie.pdf",
            #     "category_name": "全蛋白功能注释"
            # },
            "130301": {
                "main_coll": "sg_annotation_kegg",
                "name": "Pathway分类统计柱状图",
                "img": "./output/path_class.barline.pdf",
                "category_name": "全蛋白功能注释"
            },
            "130401": {
                "main_coll": "sg_annotation_kegg",
                "name": "重要通路统计图",
                "img": "./output/key_path.bar.pdf",
                "category_name": "全蛋白功能注释"
            },
            "130501": {
                "main_coll": "sg_annotation_cog",
                "name": "COG分类统计柱状图",
                "img": "./output/cog_bar.bar.pdf",
                "category_name": "全蛋白功能注释"
            },
            "130601": {
                "main_coll": "sg_annotation_pfam",
                "name": "Pfam注释柱状图",
                "img": "./output/pfam_bar.bar.pdf",
                "category_name": "全蛋白功能注释"
            },
            "130701": {
                "main_coll": "sg_annotation_subloc",
                "name": "Subloc注释柱状图",
                "img": "./output/subloc_bar.bar.pdf",
                "category_name": "全蛋白功能注释"
            },
            "140101": {
                "main_coll": "sg_diff",
                "name": "差异散点图",
                "img": "./output/*scatter.scatter.pdf",
                "category_name": "差异蛋白分析"
            },
            "140102": {
                "main_coll": "sg_diff",
                "name": "差异火山图",
                "img": "./output/*volcano.volcano.pdf",
                "category_name": "差异蛋白分析"
            },
            "150101": {
                "main_coll": "sg_proteinset_cluster",
                "name": "层级聚类热图",
                "img": "./output/cluster_group_10____cluster_group_10.heatmap.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "150201": {
                "main_coll": "sg_proteinset_cluster",
                "name": "子聚类趋势图",
                "img": "./output/cluster_group_10____sub1.line.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "151001": {
                "main_coll": "sg_proteinset_venn",
                "name": "venn图",
                "img": "./output/venn.venn.pdf",
                "category_name": "蛋白集分析"
            },
            "150301": {
                "main_coll": "sg_proteinset_go_class",
                "name": "GO分类统计柱形图",
                "img": "./output/*all_go.go_bar.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            # "150401": {
            #     "main_coll": "sg_proteinset_go_class2",
            #     "name": "GO双向直方图",
            #     "img": "./output/*up_down_go_neg.bar_neg.pdf",
            #     "category_name": "蛋白集分析",
            #     "main_dict": {"gene_set": "get_geneset_by_file_name"}
            # },
            "150501": {
                "main_coll": "sg_proteinset_go_enrich",
                "name": "GO富集分析结果图",
                "img": "./output/*go_enrich_bar.shadowbar.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            # "150502": {
            #     "main_coll": "sg_proteinset_go_enrich",
            #     "name": "GO富集分析气泡图（密集型）",
            #     "img": "./output/*go_enrich_bubble.densebubble.pdf",
            #     "category_name": "蛋白集分析",
            #     "main_dict": {"gene_set": "get_geneset_by_file_name"}
            # },
            # "150503": {
            #     "main_coll": "sg_proteinset_go_enrich",
            #     "name": "GO富集分析气泡图（分散型）",
            #     "img": "./output/*go_enrich_bubble2.scatterbubble.pdf",
            #     "category_name": "蛋白集分析",
            #     "main_dict": {"gene_set": "get_geneset_by_file_name"}
            # },
            # "150504": {
            #     "main_coll": "sg_proteinset_go_enrich",
            #     "name": "GO富集分析有向无环图",
            #     "img": "../output/5_Proteinset/04_PsetEnrich/01_EnrichGO/*/go_lineage.pdf",
            #     "category_name": "蛋白集分析",
            #     "main_dict": {"gene_set": "get_geneset_by_file_name"}
            # },
            "150701": {
                "main_coll": "sg_proteinset_kegg_class",
                "name": "Pathway分类统计柱状图",
                "img": "./output/*path_all.barline.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "150801": {
                "main_coll": "sg_proteinset_kegg_enrich",
                "name": "KEGG富集分析柱状图",
                "img": "./output/*kegg_enrich_bar.shadowbar.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            # "150802": {
            #     "main_coll": "sg_proteinset_kegg_enrich",
            #     "name": "KEGG富集分析气泡图（密集型）",
            #     "img": "./output/*kegg_enrich_bubble.densebubble.pdf",
            #     "category_name": "蛋白集分析",
            #     "main_dict": {"gene_set": "get_geneset_by_file_name"}
            # },
            # "150803": {
            #     "main_coll": "sg_proteinset_kegg_enrich",
            #     "name": "KEGG富集分析气泡图（分散型）",
            #     "img": "./output/*kegg_enrich_bubble2.scatterbubble.pdf",
            #     "category_name": "蛋白集分析",
            #     "main_dict": {"gene_set": "get_geneset_by_file_name"}
            # },
            "150901": {
                "main_coll": "sg_proteinset_ipath",
                "name": "iPath代谢通路图",
                "img": "../output/5_Proteinset/06_PsetIpath/*/Metabolism.svg",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "151101": {
                "main_coll": "sg_proteinset_circ",
                "name": "蛋白富集弦图",
                # "img": "./output/*chord.circ.pdf",
                "img": "./output/*kegg_chord.circ.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "152201": {
                "main_coll": "sg_proteinset_string_picture",
                "name": "蛋白质互作网络图(静态图)",
                "img": "../output/5_Proteinset/05_PsetStringPic/*/*.network.svg",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            # "151401": {
            #     "main_coll": "sg_proteinset_ppi",
            #     "name": "网络中心系数分布图",
            #     "img": "./output/*ppi.centrality.line.showCurve.pdf",
            #     "category_name": "蛋白集分析",
            #     "main_dict": {"gene_set": "get_geneset_by_file_name"}
            # },
            # "151501": {
            #     "main_coll": "sg_proteinset_ppi",
            #     "name": "网络节点度分布图",
            #     "img": "./output/*ppi.degree.line.showCurve.pdf",
            #     "category_name": "蛋白集分析",
            #     "main_dict": {"gene_set": "get_geneset_by_file_name"}
            # },
            "151701": {
                "main_coll": "sg_proteinset_cog_class",
                "name": "COG注释统计图",
                "img": "./output/*all_cog.bar.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "151801": {
                "main_coll": "sg_proteinset_pfam",
                "name": "Pfam注释统计柱状图",
                "img": "./output/*_all_pfam.bar.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
            },
            "152001": {
                "main_coll": "sg_proteinset_subloc",
                "name": "亚细胞定位图",
                "img": "./output/*_all_subloc.bar.pdf",
                "category_name": "蛋白集分析",
                "main_dict": {"gene_set": "get_geneset_by_file_name"}
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
            print("符合当前id：{}的文件有{}".format(str(img_id),str(files)))
            if len(files) > 1:
                if img_id == '152201':
                    pdf_choose = sorted(files)[0]
                elif img_dict["name"] not in ["GO富集分析有向无环图","iPath代谢通路图"]:
                    pdf_choose = os.path.join(os.path.split(files[0])[0], sorted([os.path.split(i)[1] for i in files])[0])
                else:
                    pdf_choose = os.path.join(os.path.split(os.path.split(files[0])[0])[0], sorted([os.path.basename(os.path.split(i)[0]) for i in files])[0], os.path.split(files[0])[1])
            elif len(files) == 1:
                pdf_choose = files[0]
            else:
                print("该分析未找到结果图片 {}".format(img_dict))
                continue
            print("最终选中的文件是{}".format(str(pdf_choose)))

            png_path = self.pdf2png(pdf_choose)
            print img_dict
            self.report_config_out[img_id] = {
                "img": png_path,
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

        # add dMass.png by hand, generated after chart
        self.report_config_out['110302'] = {
            "img": 'dMass.png',
            "pdf_choose": self.report_config['110302']['img'],
            "name": self.report_config['110302']['name'],
            "categroy_name": self.report_config['110302']["category_name"],
            "main_coll": self.report_config['110302']["main_coll"]
        }

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
        if str(pdf_choose).endswith('go_lineage.pdf'):
            png_path = os.path.splitext(os.path.basename(pdf_choose))[0]+'.png'
            command = "gs  -sDEVICE=pngalpha -o png/{} -r60 {}".format(png_path, pdf_choose)
            os.system(command)
            return png_path
        elif str(pdf_choose).endswith('pdf'):
            png_path = os.path.splitext(os.path.basename(pdf_choose))[0]
            command = "gs -o png/{}.jpg -sDEVICE=jpeg  -r256 {} && convert png/{}.jpg png/{}.png".format(png_path, pdf_choose, png_path, png_path)
            os.system(command)
            return png_path+'.png'
        elif str(pdf_choose).endswith('svg'):
            png_path = os.path.splitext(os.path.basename(pdf_choose))[0]+'.png'
            cairosvg.svg2png(url=pdf_choose, write_to=os.path.join("./png/",png_path))
            return png_path
        else:
            new_file = os.path.join("./png/",os.path.basename(pdf_choose))
            os.link(pdf_choose, new_file)
            return os.path.basename(new_file)

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
