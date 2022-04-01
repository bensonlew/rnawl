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

    def set_map_file(self, map_file):
        '''
        设置原始html路径
        '''
        self.map_file = map_file

    def set_out_file(self, out_file):
        '''
        设置原始html路径
        '''
        self.out_file = out_file

    def set_map_id(self, map_id):
        '''
        设置原始html路径
        '''
        self.map_id = map_id

    def set_title_dict(self, title_dict):
        '''
        设置KO编号与描述信息文件
        '''
        self.title_dict = title_dict

    def set_gene_dict(self, gene_dict):
        '''
        设置基因编号与描述信息文件
        '''
        self.gene_dict = gene_dict

    def change_ko_title_by_dict(self):
        '''
        重置html文件中的标签属性
        '''
        # print "+++" + self.map_file
        # print "+++" + self.map_id
        # print  self.title_dict

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
        # print self.title_dict
        for area in areas:
            area_dict=dict(area.items())
            if area_dict['shape'] == 'rect':
                if "/entry" in area_dict['href']:
                    kos = area_dict['href'].split('/')[-1].split('+')
                else:
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
        '''
        重置html文件中基因名称描述
        '''
        html = lxml.html.parse(self.map_file)
        root = html.getroot()
        body = root.find('body')

        # 修改图片链接
        img = body.find('img')
        print "+++++\n"
        # print self.map_id
        print self.gene_dict
        # try:
        #    img.set('src', self.map_id)
        # except:
        #    print "{} has no image src".format(self.map_id)

        # 修改鼠标点击矩形区域链接
        maps = body.find('map')
        areas = maps.findall('area')

        for area in areas:
            area_dict=dict(area.items())
            if area_dict['shape'] == 'rect':
                if "/entry" in area_dict['href']:
                    kos = area_dict['href'].split('/')[-1].split('+')
                else:
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

if __name__ == "__main__":
    KeggHtml = KeggHtml()
    # KeggHtml.run(map_file, out_file, png_id, title_dict)
