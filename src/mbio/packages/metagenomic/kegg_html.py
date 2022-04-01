# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET
import lxml.html

class KeggHtml(object):
    def __init__(self):
        self.map_file = ""
        self.out_file = ""
        self.map_id = ""
        self.title_dict = ""
        self.gene_dict = ""
        self.color_bg = ["#FFFF00", "#00CD00"]
        self.color_fg = ["#0000CD", "#FF0000"]
        self.c_dict = {}
        self.ko_list = []
        self.geneset_ko = []
        self.shape_list = ["rect"]
        self.all_shape = ["rect", "circle"]


    def set_map_file(self, map_file):
        """
        设置原始html路径
        """
        self.map_file = map_file

    def set_ko_list(self, ko):
        """
        设置已知基因与新基因列表
        """
        self.ko_list = ko

    def set_out_file(self, out_file):
        """
        设置原始html路径
        """
        self.out_file = out_file

    def set_map_id(self, map_id):
        """
        设置原始html路径
        """
        self.map_id = map_id

    def set_gene_dict(self, gene_dict):
        """
        设置基因编号与描述信息文件
        """
        self.gene_dict = gene_dict

    def set_geneset_ko(self, geneset_ko):
        """
        设置基因集与ko列表
        """
        self.geneset_ko = geneset_ko

    def set_shape_list(self, shape_list=None):
        """
        设置需要处理的图形元素形状，rect/circle
        :param shape_list: ["rect", "circle"]
        """
        if shape_list:
            self.shape_list = shape_list

    def get_bg_color(self, kos):
        """
        获取标记背景色
        """
        if len(self.ko_list) == 1:
            return self.color_bg[0]
        elif len(self.ko_list) == 2:
            colors = []
            if set(kos).intersection(self.ko_list[0]):
                colors.append(self.color_bg[0])
            if set(kos).intersection(self.ko_list[1]):
                colors.append(self.color_bg[1])
            return ",".join(colors)
        else:
            raise Exception("ko 列表只能包含一个或两个set")

    def get_fg_color(self, kos):
        """
        获取标记边框色
        """
        if len(self.geneset_ko) == 1:
            if len(set(kos).intersection(self.geneset_ko[0])) >= 1:
                return self.color_fg[0]
            else:
                return ""
        elif len(self.geneset_ko) == 2:
            colors = []
            if set(kos).intersection(self.geneset_ko[0]):
                colors.append(self.color_fg[0])
            if set(kos).intersection(self.geneset_ko[1]):
                colors.append(self.color_fg[1])
            return ",".join(colors)
        else:
            raise Exception("geneset 列表只能包含一个或两个set")

    def get_geneset_title(self, ko_title):
        """
        根据基因集修改ko_title描述信息
        """
        titles = ko_title.split('\\n')
        titles_new = []
        for title in titles:
            if "accession: " in title:
                # 蛋白产品显示格式
                titles_ele = title.split("accession: ")[1]
                genes = titles_ele.split("|")
                genes_new = []
                for gene in genes:
                    if self.gene_dict.has_key(gene):
                        genes_new.append(gene + "(" + self.gene_dict[gene] + ")")
                    else:
                        genes_new.append(gene)
                title_new = title.split("accession: ")[0] + "accession: " + "|".join(genes_new)
            if ":" in title:
                # rna产品展示格式
                titles_ele = title.split(":")[1]
                titles_ele = titles_ele.rstrip(";")
                genes = titles_ele.split(",")
                genes_new = []
                for gene in genes:
                    if self.gene_dict.has_key(gene):
                        genes_new.append(gene + "(" + self.gene_dict[gene] + ")")
                    else:
                        genes_new.append(gene)
                title_new = title.split(":")[0] + ":" + ",".join(genes_new)
            else:
                title_new = title
            titles_new.append(title_new)
        all_title_new = "\\n".join(titles_new)
        return all_title_new

    def read_html(self):
        """
        读取html文件
        :return: areas
        """
        html = lxml.html.parse(self.map_file)
        root = html.getroot()
        body = root.find('body')

        # 修改图片链接
        img = body.find('img')
        try:
            img.set('src', self.map_id)
        except:
            print "{} has no image src".format(self.map_id)

        # 修改鼠标点击矩形区域链接
        maps = body.find('map')
        areas = maps.findall('area')
        return html,areas

    def change_ko_title_by_dict(self):
        """
        重置html文件中的标签属性
        """
        html,areas = self.read_html()
        for area in areas:
            area_dict=dict(area.items())
            if area_dict['shape'] in self.shape_list:
                kos = area_dict['href'].split('?')[-1].split('+')
                # print "***".join(kos)
                change_kos = list(set(kos) & set(self.title_dict.keys()))
                # print "change***".join(change_kos)
                if len(change_kos) > 0:
                    titles = area_dict['title'].split(',')
                    titles_new = []
                    for title in titles:
                        title_new = title + ":"
                        for ko in change_kos:
                            if ko in title:
                                title_new += self.title_dict[ko] + ";"
                        titles_new.append(title_new)

                    all_title_new = "\n".join(titles_new)
                    # print "***" + all_title_new
                    area.set('title', all_title_new)
        html.write(self.out_file)

    def change_gene_by_dict(self):
        """
        重置html文件中基因名称描述
        """
        html,areas = self.read_html()
        for area in areas:
            area_dict=dict(area.items())
            if area_dict['shape'] in self.shape_list:
                kos = area_dict['href'].split('?')[-1].split('+')
                # print "***".join(kos)
                if 1:
                    titles = area_dict['title'].split('\n')
                    titles_new = []
                    for title in titles:
                        if "accession: " in title:
                            # 蛋白产品显示格式
                            titles_ele = title.split("accession: ")[1]
                            genes = titles_ele.split("|")
                            genes_new = []
                            for gene in genes:
                                if self.gene_dict.has_key(gene):
                                    genes_new.append(gene + "(" + self.gene_dict[gene] + ")")
                                else:
                                    genes_new.append(gene)
                            title_new = title.split("accession: ")[0] + "accession: " + "|".join(genes_new)
                        if ":" in title:
                            # rna产品展示格式
                            titles_ele = title.split(":")[1]
                            titles_ele = titles_ele.rstrip(";")
                            genes = titles_ele.split(",")
                            genes_new = []
                            for gene in genes:
                                if self.gene_dict.has_key(gene):
                                    genes_new.append(gene + "(" + self.gene_dict[gene] + ")")
                                else:
                                    genes_new.append(gene)
                            title_new = title.split(":")[0] + ":" + ",".join(genes_new)
                        else:
                            title_new = title
                        titles_new.append(title_new)
                    all_title_new = "\n".join(titles_new)
                    # print "***" + all_title_new
                    area.set('title', all_title_new)
        html.write(self.out_file)

    def get_html_mark(self, reset_title=False):
        """
        获取标签信息，导入用于前端作图
        """
        html,areas = self.read_html()
        with open(self.out_file, 'w') as mark:
            for area in areas:
                area_dict=dict(area.items())
                if area_dict['shape'] in self.shape_list:
                    kos = area_dict['href'].split('?')[-1].split('+')
                    # print "***".join(kos)
                    change_kos = list(set(kos) & set(self.title_dict.keys()))
                    # print "change***".join(change_kos)
                    if len(change_kos) > 0:
                        titles = area_dict['title'].split(',')
                        if reset_title:
                            titles_new = []
                            for title in titles:
                                title_new = title + ":"
                                for ko in change_kos:
                                    if ko in title:
                                        title_new += self.title_dict[ko] + ";"
                                titles_new.append(title_new)

                            all_title_new = "\n".join(titles_new)
                            all_title_mark = "\\n".join(titles_new)
                            # print "***" + all_title_new
                            area.set('title', all_title_new)
                        else:
                            all_title_mark = area_dict['title']
                        bg_color = self.get_bg_color(change_kos)
                        mark.write("\t".join([self.map_id, area_dict['shape'], bg_color, "",  area_dict['coords'],  all_title_mark, ",".join(change_kos), area_dict['href']]) + "\n")
                    else:
                        mark.write("\t".join([self.map_id, area_dict['shape'], "", "",  area_dict['coords'],  area_dict['title'], "", area_dict['href']]) + "\n")
                elif area_dict['shape'] in self.all_shape:
                    mark.write("\t".join([self.map_id, area_dict['shape'], "", "",  area_dict['coords'],  area_dict['title'], "", area_dict['href']]) + "\n")

    def get_mark_from_html(self, html):
        """
        使用html文件获取mark标记， 用于之前项目中只有html文件， 没有mark标记的情况
        """
        html,areas = self.read_html()
        with open("temp.mark", 'w') as mark:
            for area in areas:
                area_dict=dict(area.items())
                if area_dict['shape'] in self.shape_list:
                    kos = area_dict['href'].split('?')[-1].split('+')
                    # print "***".join(kos)
                    change_kos = list(set(kos) & set(self.title_dict.keys()))
                    all_title = area_dict['title'].replace("\n", "\\n")
                    bg_color = self.get_bg_color(change_kos)
                    mark.write("\t".join([self.map_id, area_dict['shape'], bg_color, "",  area_dict['coords'],  all_title, ",".join(change_kos), area_dict['href']]) + "\n")
                elif area_dict['shape'] in self.all_shape:
                    mark.write("\t".join([self.map_id, area_dict['shape'], "", "",  area_dict['coords'],  area_dict['title'], "", area_dict['href']]) + "\n")
        return "temp.mark"

    def get_mark_from_html2(self, html):
        """
        使用html文件获取mark标记， 用于在删除了mark标记的情况下恢复文件
        """
        html, areas = self.read_html()
        with open("temp.mark", 'w') as mark:
            for area in areas:
                area_dict=dict(area.items())
                if area_dict['shape'] in self.shape_list:
                    kos = area_dict['href'].split('?')[-1].split('+')
                    # print "***".join(kos)
                    # change_kos = list(set(kos) & set(self.title_dict.keys()))
                    change_kos = []
                    all_title = area_dict['title'].replace("\n", "\\n")

                    bg_color = ""
                    if ":" in all_title:
                        if "MSTRG" in all_title or "XLOC" in all_title or "TCONS" in all_title:
                            bg_color = self.color_bg[0] + "," + self.color_bg[1]
                        else:
                            bg_color = self.color_bg[0]

                        titles = all_title.split("\\n")
                        for title in titles:
                            if ";" in title:
                                change_kos.append(title.strip().split()[0])
                    else:
                        pass
                    # bg_color = self.get_bg_color(change_kos)
                    mark.write("\t".join([self.map_id, area_dict['shape'], bg_color, "",  area_dict['coords'],  all_title, ",".join(change_kos), area_dict['href']]) + "\n")
                elif area_dict['shape'] in self.all_shape:
                    mark.write("\t".join([self.map_id, area_dict['shape'], "", "",  area_dict['coords'],  area_dict['title'], "", area_dict['href']]) + "\n")
        return "temp.mark"

    def get_gene_set_html_mark(self, html_mark):
        """
        获取标签信息，导入用于前端作图
        """
        with open(html_mark, 'r') as mark_in, open(self.out_file, 'w') as mark_out:
            for line in mark_in.readlines():
                cols = line.strip().split("\t")
                kos = cols[6].split(",")
                if cols[1] == 'rect':
                    cols[3] = self.get_fg_color(kos)
                    cols[5] = self.get_geneset_title(cols[5])
                    mark_out.write("\t".join(cols) + "\n")
                elif cols[1] == "circle":
                    mark_out.write("\t".join(cols) + "\n")

    def get_mg_mark(self):
        """
        根据self.ko_list 与 map_file获取唯一的蛋白编号，或酶，酶优先，多个取第一
        ko_list中除蛋白，可能包含酶
        :return:
        """
        html,areas = self.read_html()
        kos_list = []
        enzyme_list = []
        with open(self.out_file, 'w') as mark:
            for area in areas:
                tmp_ko = ""
                area_dict = dict(area.items())
                if area_dict['shape'] in self.shape_list:
                    kos = area_dict['href'].split('?')[-1].split('+')
                    change_kos = self.geneset_ko[0] & set(kos)
                    print kos
                    beixuan_kos = []
                    for one in kos:
                        if one.startswith("K") and one in change_kos:
                            beixuan_kos.append(one)  # 保留蛋白编号和href中的顺序
                        elif one in change_kos:
                            enzyme_list.append(one)  # 优先使用酶编号
                            tmp_ko = one
                            if len(beixuan_kos) != 0:
                                break
                    if len(enzyme_list) == 0 and len(change_kos) != 0:  # 没有对应酶，用对应的第一个蛋白代替
                        kos_list.append(beixuan_kos[0])
                        tmp_ko = beixuan_kos[0]
                    # fg_color = self.get_fg_color(kos)
                    bg_color = self.color_bg[0]  if tmp_ko else tmp_ko # 对于kegg注释，设置此处做背景色，如果是代谢通路差异分析，此处设置的颜色需忽略
                    mark.write("\t".join([self.map_id, area_dict['shape'], bg_color, "",  area_dict['coords'],  area_dict['title'], tmp_ko, area_dict['href']]) + "\n")
                elif area_dict['shape'] in self.all_shape:
                    mark.write("\t".join([self.map_id, area_dict['shape'], "", "",  area_dict['coords'],  area_dict['title'], "", area_dict['href']]) + "\n")
        return set(kos_list),set(enzyme_list)

    def run(self, map_file, out_file, png_id, title_dict):
        self.set_map_file(map_file)
        self.set_out_file(out_file)
        self.set_map_id(png_id)
        self.set_title_dict(title_dict)
        self.change_ko_title_by_dict()

    def run_gene_set(self, map_file, out_file, gene_dict):
        self.set_map_file(map_file)
        self.set_out_file(out_file)
        self.set_gene_dict(gene_dict)
        self.change_gene_by_dict()

    def run_html_mark(self, map_file, out_file, png_id, title_dict, ko):
        self.set_map_file(map_file)
        self.set_out_file(out_file)
        self.set_map_id(png_id)
        self.set_title_dict(title_dict)
        self.set_ko_list(ko)
        self.get_html_mark(reset_title=False)

    def run_gene_set_mark(self, map_file, out_file, gene_dict, html_mark, geneset_ko):
        if os.path.exists(html_mark):
            self.set_map_file(map_file)
            
            self.set_out_file(out_file)
            self.set_gene_dict(gene_dict)
            # print geneset_ko
            self.set_geneset_ko(geneset_ko)
            self.get_gene_set_html_mark(html_mark)
        else:
            return 1

    def run_mg_mark(self, map_file, out_file, png_id, geneset_ko):
        """
        宏基因组kegg通路差异分析, 对于kegg注释，整个注释相当于一个geneset_ko
        """
        self.set_map_file(map_file)
        self.set_out_file(out_file)
        self.set_map_id(png_id)
        self.set_geneset_ko(geneset_ko)
        kos_list,enzyme_list = self.get_mg_mark()
        return kos_list,enzyme_list

    def run_gene_set_mark_from_html(self, map_id, map_file, out_file, gene_dict, html, geneset_ko):
        # html文件里缺少信息暂时无法兼容
        if os.path.exists(html):
            self.set_map_id(map_id)
            self.set_map_file(map_file)
            self.set_out_file(out_file)
            self.set_gene_dict(gene_dict)
            self.set_geneset_ko(geneset_ko)
            html_mark = self.get_mark_from_html2(html)
            self.get_gene_set_html_mark(html_mark)
        else:
            return 1


if __name__ == "__main__":
    KeggHtml = KeggHtml()
    # KeggHtml.run(map_file, out_file, png_id, title_dict)
