# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET

class Ipath(object):
    def __init__(self):
        project_type = "ref_rna"
        self.db_dir = ""
        self.maps = ['Metabolic_pathways', 'Regulatory_pathways', 'Biosynthesis_of_secondary_metabolities']
        self.K = dict()
        self.legend = []

    def set_db(self, db_dir):
        '''
        设置ipath数据库路径，包括下载的svg文件和map图形KO对应关系
        '''
        self.db_dir = db_dir

    def get_K_color_width(self, color_width_file):
        '''
        获得KO和线条长度颜色信息
        '''
        with open(color_width_file, 'r') as f:
            for line in f.readlines():
                k, color, width = line.strip().split('\t')
                width=width.replace('W','')
                if self.K.has_key(k):
                    pass
                else:
                    self.K.update({k:{"color":color, "width":width}})
    def set_legend(self, sets):
        '''
        添加基因集标签
        '''
        self.legend.extend(sets)

    def map_ko_to_svg(self, K, id_file):
        '''
        将KO的颜色长度信息映射至图形svg
        '''
        svg_label = dict()
        with open(id_file, 'r') as f:
            for line in f.readlines():
                label, first_K, all_K = line.strip().split('\t')
                all_ele = [x.strip().split(" ")[0] for x in all_K.split(",")]
                for ele in all_ele:
                    if K.has_key(ele):
                        if svg_label.has_key(label):
                            if K[ele]["color"] == "#0000ff":
                                svg_label.update({label:{"color":K[ele]["color"],"width":K[ele]["width"]}})
                            # pass
                        else:
                            svg_label.update({label:{"color":K[ele]["color"],"width":K[ele]["width"]}})
                    else:
                        pass
        return svg_label

    def set_svg_legend(self, element):
        '''
        修改图例
        '''
        leg_col = element.find('{http://www.w3.org/2000/svg}path')
        text=element.find('{http://www.w3.org/2000/svg}text')
        text.text = self.legend[0]
        text_x = int(text.attrib['x'])
        text_font_size = int(text.attrib['style'].split(";")[0].strip()[-4:-2])
        if len(self.legend) == 2:
            keg_col2= leg_col.copy()
            text2 = text.copy()
            keg_col3= leg_col.copy()
            text3 = text.copy()

            text2_x = text_x  + int(text_font_size * 0.6 * len(self.legend[0])) + 40
            keg_col2.attrib['d'] = 'M{},30L{},30'.format(text2_x , text2_x+100)
            text2.attrib['x'] = str(text2_x + 150)
            text2.attrib['fill'] = "#00ff00"
            keg_col2.attrib['stroke'] = "#00ff00"

            text3_x = text2_x + 150 + int(text_font_size * 0.6 * len(self.legend[1])) + 40
            keg_col3.attrib['d'] = 'M{},30L{},30'.format(text3_x , text3_x+100)
            text3.attrib['x'] = str(text3_x + 150)

            keg_col3.attrib['stroke'] = "#0000ff"
            text3.attrib['fill'] = "#0000ff"

            text2.text = self.legend[1]
            # print "legend 2 is {}".format(self.legend[1])
            text3.text = 'both'

            element.insert(2, keg_col2)
            element.insert(3, text2)
            element.insert(4, keg_col3)
            element.insert(5, text3)

    def map_file(self):
        '''
        修改svg图形改变其中的线条颜色和宽度
        '''
        for path_map in self.maps:
            map_file = os.path.join(self.db_dir, path_map + ".raw.svg")
            id_file = os.path.join(self.db_dir, path_map + ".map")

            svg_lable = self.map_ko_to_svg(self.K, id_file)
            out_file = path_map + ".svg"
            ET.register_namespace("","http://www.w3.org/2000/svg")
            xml = ET.parse(map_file)
            root = xml.getroot()
            layer_a = root.findall("{http://www.w3.org/2000/svg}g")
            for one_of_a in layer_a:
                if one_of_a.attrib['id'] == 'title':
                    continue
                if one_of_a.attrib['id'] == 'legend':
                    self.set_svg_legend(one_of_a)
                    continue
                layer_b = one_of_a.findall("{http://www.w3.org/2000/svg}g")
                for one_of_b in layer_b:
                    # 修改标签

                    # print self.legend
                    if one_of_b.attrib['id'] == 'title':
                        continue
                    if one_of_b.attrib['id'] == 'legend':
                        self.set_svg_legend(one_of_b)
                        continue

                    layer_c = one_of_b.findall("{http://www.w3.org/2000/svg}g")
                    for one_of_c in layer_c:
                        m_id = "no_id"
                        if one_of_c.attrib.has_key('id'):
                            m_id = one_of_c.attrib['id']
                        else:
                            pass
                        layer_d = one_of_c.findall("{http://www.w3.org/2000/svg}path")
                        layer_d2 = one_of_c.findall("{http://www.w3.org/2000/svg}ellipse")
                        for one_of_d in layer_d + layer_d2:
                            # print m_id,
                            attr_list = ['style', 'stroke-width' , 'stroke','fill']
                            if svg_lable.has_key(m_id) and one_of_d.attrib.has_key("stroke-width"):
                                one_of_d.attrib["stroke-width"] = svg_lable[m_id]["width"]
                                one_of_d.attrib["style"] =  ' opacity: 0.8; '
                            else:
                                pass

                            if svg_lable.has_key(m_id) and one_of_d.attrib.has_key("stroke"):
                                one_of_d.attrib["stroke"] = svg_lable[m_id]["color"]
                                if svg_lable[m_id]["color"] == "#0000ff":
                                    one_of_d.attrib["style"] =  ' opacity: 0.8; '
                                elif len(self.legend) == 2:
                                    one_of_d.attrib["style"] =  ' opacity: 0.4; '
                            else:
                                pass

            xml.write(out_file)

if __name__ == "__main__":
    Ipath1 = Ipath()
    Ipath1.set_db(sys.argv[1])
    Ipath1.get_K_color_width(sys.argv[2])
    Ipath1.map_file()
