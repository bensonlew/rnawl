# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import re
import os
import sys
import xml.etree.ElementTree as ET
import pandas as pd


class Ipath(object):
    def __init__(self):
        self.maps = ['Metabolic_pathways', 'Regulatory_pathways', 'Biosynthesis_of_secondary_metabolities']
        self.K = dict()
        self.db_dir = ""
        self.legend = []
        self.share_color = ""
        self.svg_label = dict()  # map_ko_to_svg时使用
        self.color_map = dict()  # 各个legend对应的颜色

    def set_db(self, db_dir):
        """
        设置ipath数据库路径，包括下载的svg文件和map图形KO对应关系
        :param db_dir:
        :return:
        """
        self.db_dir = db_dir

    def set_legend(self, sets):
        """
        添加基因集标签
        :param sets:
        :return:
        """
        self.legend.extend(sets)

    def set_share_color(self, color):
        """
        所有分组共有是的颜色
        :return:
        """
        self.share_color = color

    def set_color_map(self, color_map):
        self.color_map = color_map
        self.share_color = color_map.get("Share", self.share_color)

    def get_K_color_width(self, color_width_file):
        """
        从ipath_input.xls获取KO对应线条长度和颜色的信息
        :param color_width_file:
        :return:
        """
        cwf = pd.read_table(color_width_file, names=["color", "width"], index_col=0)
        cwf.dropna(axis=0, how="any", inplace=True)
        cwf.replace({"W15": "15"}, inplace=True)
        self.K = cwf.to_dict("index")

    def parse_svg(self, df):
        label_list = df['label'].split(",")
        if self.K.has_key(df['pro']):
            pro_info = self.K[df['pro']]
            for label in label_list:
                # if pro_info['color'] != self.share_color:
                #    continue
                if self.svg_label.has_key(label) and self.svg_label[label]["color"] == self.share_color:
                    continue  # 防止share颜色被覆盖
                self.svg_label.update({label: {"color": pro_info["color"], "width": pro_info["width"]}})


    def map_ko_to_svg(self, id_file):
        """
        将KO颜色，长度信息映射至图形svg
        :param id_file:
        :return:
        """
        self.svg_label = dict()  # 每次运行需要清空
        data = pd.read_table(id_file, names=['pro', 'label'], index_col=None)
        data.apply(self.parse_svg, axis=1)

    def set_svg_legend(self, element, root):
        """
        修改图例
        :param element:
        :param root: 需要获取width,height,并修改height
        :return:
        """
        width = int(root.attrib["width"])
        height = int(root.attrib["height"])
        leg_col = element.find('{http://www.w3.org/2000/svg}path')
        text=element.find('{http://www.w3.org/2000/svg}text')
        resort_legend = self.legend[1:] + self.legend[:1]  # 将Share组放到最后
        # text_x = int(text.attrib['x'])
        text_x = 30
        text_y = height + 40
        text.text = resort_legend[0]
        text.attrib['x'] = str(text_x + 150)
        text.attrib['y'] = str(text_y + 15)
        text.attrib['fill'] = self.color_map[resort_legend[0]]
        leg_col.attrib['stroke'] = self.color_map[resort_legend[0]]
        leg_col.attrib['d'] = 'M{0},{1}L{2},{1}'.format(text_x, text_y, text_x+100)
        text_font_size = int(text.attrib['style'].split(";")[0].strip()[-4:-2])
        text_x = text_x + 150 + int(text_font_size * 0.6 * len(resort_legend[0])) + 40
        for index,group in enumerate(resort_legend):
            if index == 0:
                continue
            keg_col = leg_col.copy()
            keg_text = text.copy()
            text_x_border = text_x + 150 + int(text_font_size * 0.6 * len(group)) + 40
            if text_x_border > width:
                text_x = 30
                text_y += 70
            keg_col.attrib['d'] = 'M{0},{1}L{2},{1}'.format(text_x , text_y, text_x+100)
            keg_text.attrib['x'] = str(text_x + 150)
            keg_text.attrib['y'] = str(text_y + 15)
            keg_text.attrib['fill'] = self.color_map[group]
            keg_col.attrib['stroke'] = self.color_map[group]
            text_x = text_x + 150 + int(text_font_size * 0.6 * len(group)) + 40
            keg_text.text = group
            element.insert(2*(index+1), keg_col)
            element.insert(2*(index+1)+1, keg_text)
        root.attrib["height"] = str(text_y + 70)
        box_coord = root.attrib["viewBox"].split(" ")
        box_coord[3] = str(text_y + 90)
        root.attrib["viewBox"] = " ".join(box_coord)

    def map_file(self):
        for path_map in self.maps:
            map_file = os.path.join(self.db_dir, path_map + ".raw.svg")
            id_file = os.path.join(self.db_dir, path_map + ".metag.map")
            self.map_ko_to_svg(id_file)
            out_file = path_map + ".svg"
            xml = ET.parse(map_file)
            root = xml.getroot()
            layer_a = root.findall("{http://www.w3.org/2000/svg}g")
            for one_of_a in layer_a:
                if one_of_a.attrib['id'] == 'title':
                    continue
                if one_of_a.attrib['id'] == 'legend':
                    self.set_svg_legend(one_of_a, root)
                    continue
                layer_b = one_of_a.findall("{http://www.w3.org/2000/svg}g")
                for one_of_b in layer_b:
                    if one_of_b.attrib['id'] == 'title':
                        continue
                    if one_of_b.attrib['id'] == 'legend':
                        self.set_svg_legend(one_of_b, root)
                        continue
                    layer_c = one_of_b.findall("{http://www.w3.org/2000/svg}g")
                    for one_of_c in layer_c:
                        m_id = "no_id"
                        if one_of_c.attrib.has_key('id'):
                            m_id = one_of_c.attrib['id']
                        layer_d = one_of_c.findall("{http://www.w3.org/2000/svg}path")
                        layer_d2 = one_of_c.findall("{http://www.w3.org/2000/svg}ellipse")
                        for one_of_d in layer_d + layer_d2:
                            if self.svg_label.has_key(m_id) and one_of_d.attrib.has_key("stroke-width"):
                                one_of_d.attrib["stroke-width"] = self.svg_label[m_id]["width"]
                                one_of_d.attrib["style"] =  ' opacity: 1; '
                            if self.svg_label.has_key(m_id) and one_of_d.attrib.has_key("stroke"):
                                one_of_d.attrib['stroke'] = self.svg_label[m_id]["color"]
                                if self.svg_label[m_id]['color'] == self.share_color:
                                    one_of_d.attrib["style"] = ' opacity: 1; '
                                elif len(self.legend) > 2:
                                    one_of_d.attrib["style"] =  ' opacity: 0.5; '
            xml.write(out_file)