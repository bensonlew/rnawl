# -*- coding: utf-8 -*-


import re
import os
import sys
import xml.etree.ElementTree as ET
import subprocess

class Ipath(object):
    def __init__(self):
        project_type = "metabolome"
        self.db_dir = ""
        self.maps = ['Metabolism', 'Secondary_metabolites', 'Antibiotics', 'Microbial_metabolism']

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
         K：有compound的颜色和width，id_file：有compound和svg id的对应信息。据此配置svg id的颜色和width
        '''
        svg_label = dict()
        with open(id_file, 'r') as f:
            for line in f.readlines():
                spline = line.strip().split("\t")
                if len(spline) == 2:
                    compound, label = line.strip().split("\t")
                else:
                    continue
                label_list = label.split(",")
                if K.has_key(compound) and label_list[0] != "-":
                    for label in label_list:
                        if svg_label.has_key(label):
                            if K[compound]["color"] == "#0000ff":
                                svg_label.update({label:{"color":K[compound]["color"],"width":K[compound]["width"]}})
                        else:
                            svg_label.update({label: {"color":K[compound]["color"],"width":K[compound]["width"]}})

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
            all_len = len(self.legend[0]) + len(self.legend[1]) + 4
            text_font_size_new = int(2000/all_len)
            text_font_size_new = min(50, text_font_size_new)

            text.attrib['style'] = text.attrib['style'].replace(str(text_font_size) + "px", str(text_font_size_new) + "px")
            
            keg_col2= leg_col.copy()
            text2 = text.copy()
            keg_col3= leg_col.copy()
            text3 = text.copy()

            text2_x = text_x  + int(text_font_size_new * 0.6 * len(self.legend[0])) + 40
            keg_col2.attrib['d'] = 'M{},75L{},75'.format(text2_x , text2_x+100)
            text2.attrib['x'] = str(text2_x + 150)
            text2.attrib['fill'] = "#00ff00"
            keg_col2.attrib['stroke'] = "#00ff00"

            text3_x = text2_x + 150 + int(text_font_size_new * 0.6 * len(self.legend[1])) + 40
            keg_col3.attrib['d'] = 'M{},75L{},75'.format(text3_x , text3_x+100)
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

    def map_file(self,image_magick):
        '''
        修改svg图形改变其中的线条颜色和宽度
        '''
        for path_map in self.maps:
            map_file = os.path.join(self.db_dir, path_map + ".raw.svg")
            id_file = os.path.join(self.db_dir, path_map+'.compound.map')  # compound id
            k_file = os.path.join(self.db_dir, path_map+'.k.map')  # KO id

            svg_lable = self.map_ko_to_svg(self.K, id_file)
            svg_lable2 = self.map_ko_to_svg(self.K, k_file)
            svg_lable.update(svg_lable2)

            out_file = path_map + ".svg"
            out_file_png = path_map + ".png"
            xml = ET.parse(map_file)
            root = xml.getroot()
            # layer_a = root.findall("{http://www.w3.org/2000/svg}g")
            layer_a = root.findall("{http://www.w3.org/2000/svg}g")[0].findall("{http://www.w3.org/2000/svg}g")[0].findall("{http://www.w3.org/2000/svg}g")
            # print layer_a

            for one_of_a in layer_a[0].findall("{http://www.w3.org/2000/svg}g"):
                if one_of_a.attrib['id'] == 'title':
                    continue
                if one_of_a.attrib['id'] == 'legend':
                    self.set_svg_legend(one_of_a)
                    continue

            layer_b = layer_a[1].findall("{http://www.w3.org/2000/svg}g")
            for one_of_b in layer_b:
                # 修改标签
                '''
                # print self.legend
                if one_of_b.attrib['id'] == 'title':
                    continue
                if one_of_b.attrib['id'] == 'legend':
                    self.set_svg_legend(one_of_b)
                    continue

                layer_c = one_of_b.findall("{http://www.w3.org/2000/svg}g")
                for one_of_c in layer_c:
                '''
                # print "***", one_of_b
                # m_id = "no_id"
                # if one_of_b.attrib.has_key('id'):
                #     m_id = one_of_b.attrib['id']
                # else:
                #     pass
                #layer_d = one_of_b.findall("{http://www.w3.org/2000/svg}path")
                layer_d2 = one_of_b.findall("{http://www.w3.org/2000/svg}ellipse")

                #for one_of_d in layer_d + layer_d2:
                for one_of_d in layer_d2:
                    # print m_id,
                    attr_list = ['style', 'stroke-width' , 'stroke', 'fill']
                    if 'id' in one_of_d.attrib:
                        m_id = one_of_d.attrib['id']

                    if svg_lable.has_key(m_id) and one_of_d.attrib.has_key("stroke"):

                        # print "****", m_id
                        one_of_d.attrib["stroke"] = svg_lable[m_id]["color"]
                        style = one_of_d.attrib['style']
                        style_dict = {k.strip(): v.strip() for k, v in (item.strip().split(':') for item in style.strip(';').split(';'))}
                        style_dict['stroke-width'] = svg_lable[m_id]["width"]
                        if svg_lable[m_id]["color"] == "#0000ff":
                            style_dict['opacity'] = 0.8
                        elif len(self.legend) == 2:
                            style_dict['opacity'] = 0.4
                        one_of_d.attrib['style'] = "; ".join(["{}: {}".format(k, v) for k,v in style_dict.items()])
                        print  one_of_d.attrib

                    if  svg_lable.has_key(m_id) and one_of_d.attrib.has_key("fill"):
                        one_of_d.attrib["fill"] = svg_lable[m_id]["color"]
                        if one_of_d.attrib.has_key("rx") and one_of_d.attrib.has_key("ry"):
                            one_of_d.attrib["rx"] = '10'
                            one_of_d.attrib["ry"] = '10'

                    # else:
                    #     if one_of_d in layer_d:
                    #         one_of_d.attrib["stroke"] = "#999999"
                    #     else:
                    #         one_of_d.attrib["fill"] = "#999999"
                    #     pass

            xml.write(out_file)
            cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + out_file + ' ' + out_file_png
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError:
                print "svg trans to png faied"

if __name__ == "__main__":
    Ipath1 = Ipath()
    Ipath1.set_db(sys.argv[1])
    Ipath1.set_legend(set(["test"]))
    Ipath1.get_K_color_width(sys.argv[2])
    image_magick = sys.argv[3]
    Ipath1.map_file(image_magick)
