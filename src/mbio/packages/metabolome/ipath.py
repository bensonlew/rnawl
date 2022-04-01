# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import re
import os
import sys,subprocess
import xml.etree.ElementTree as ET
from mbio.packages.itraq_and_tmt.ipath import Ipath

class Ipath(Ipath):
    def __init__(self):
        super(Ipath, self).__init__()
        project_type = "metabolome"
        self.maps = ['Metabolic_pathways', 'Regulatory_pathways', 'Biosynthesis_of_secondary_metabolities']

    def map_ko_to_svg(self, K, id_file):
        '''
        将KO的颜色长度信息映射至图形svg
        '''
        svg_label = dict()
        with open(id_file, 'r') as f:
            for line in f.readlines():
                # label, first_K, all_K = line.strip().split('\t')
                compound, label = line.strip().split("\t")
                # all_ele = [x.strip().split(" ")[0] for x in all_K.split(",")]
                label_list = label.split(",")
                if K.has_key(compound) and label_list[0] != "-":
                    for label in label_list:
                        if svg_label.has_key(label):
                            if K[compound]["color"] == "#0000ff":
                                svg_label.update({label:{"color":K[compound]["color"],"width":K[compound]["width"]}})
                        else:
                            svg_label.update({label: {"color":K[compound]["color"],"width":K[compound]["width"]}})
        return svg_label

    def map_file(self, image_magick):
        '''
        修改svg图形改变其中的线条颜色和宽度
        '''
        for path_map in self.maps:
            map_file = os.path.join(self.db_dir, path_map + ".raw.svg")
            id_file = os.path.join(self.db_dir, path_map + ".compound.map")
            svg_lable = self.map_ko_to_svg(self.K, id_file)
            out_file = path_map + ".svg"
            out_file_png = path_map + ".png"
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
            cmd = image_magick + ' -flatten -quality 100 -density 130 -background white ' + out_file + ' ' + out_file_png
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError:
                print "svg trans to png faied"

if __name__ == "__main__":
    Ipath1 = Ipath()
    Ipath1.set_db(sys.argv[1])
    Ipath1.get_K_color_width(sys.argv[2])
    Ipath1.map_file()
